/* recovered: for all sites not present in a counts file, add those sites as
 * non-covered and with the appropriate context.
 *
 * Copyright (C) 2023 Andrew D. Smith
 *
 * Authors: Andrew D. Smith
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 */

#include "MSite.hpp"
#include "bsutils.hpp"
#include "counts_header.hpp"

#include <bamxx.hpp>

// from smithlab_cpp
#include "OptionParser.hpp"
#include "smithlab_os.hpp"
#include "smithlab_utils.hpp"

#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

static void
verify_chrom_orders(
  const bool verbose, const std::uint32_t n_threads,
  const std::string &filename,
  const std::unordered_map<std::string, std::int32_t> &chroms_order) {
  bamxx::bgzf_file in(filename, "r");
  if (!in)
    throw std::runtime_error("bad file: " + filename);

  bamxx::bam_tpool tp(n_threads);
  // set the threads for the input file decompression
  if (n_threads > 1 && in.is_bgzf())
    tp.set_io(in);

  std::unordered_set<std::int32_t> chroms_seen;
  std::string line;
  std::string prev_chrom;

  std::int32_t prev_idx = -1;

  while (getline(in, line)) {
    if (is_counts_header_line(line))
      continue;
    line.resize(line.find_first_of(" \t"));
    if (line != prev_chrom) {
      if (verbose)
        std::cerr << "verifying: " << line << "\n";

      const auto idx_itr = chroms_order.find(line);
      if (idx_itr == cend(chroms_order))
        throw std::runtime_error("chrom not found genome file: " + line);
      const auto idx = idx_itr->second;

      if (chroms_seen.find(idx) != end(chroms_seen))
        throw std::runtime_error("chroms out of order in: " + filename);
      chroms_seen.insert(idx);

      if (idx < prev_idx)
        throw std::runtime_error("inconsistent chromosome order at: " + line);

      prev_idx = idx;
      std::swap(line, prev_chrom);
    }
  }
  if (verbose)
    std::cerr << "chrom orders are consistent\n";
}

struct quick_buf : public std::ostringstream,
                   public std::basic_stringbuf<char> {
  // ADS: By user ecatmur on SO; very fast. Seems to work...
  quick_buf() {
    // ...but this seems to depend on data layout
    static_cast<std::basic_ios<char> &>(*this).rdbuf(this);
  }
  void
  clear() {
    // reset buffer pointers (member functions)
    setp(pbase(), pbase());
  }
  char const *
  c_str() {
    /* between c_str and insertion make sure to clear() */
    *pptr() = '\0';
    return pbase();
  }
};

/* The three functions below here should probably be moved into
   bsutils.hpp. I am not sure if the DDG function is needed, but it
   seems like if one considers strand, and the CHH is not symmetric,
   then one needs this. Also, Qiang should be consulted on this
   because he spent much time thinking about it in the context of
   plants. */
static inline bool
is_chh(const std::string &s, std::size_t i) {
  return (i < (s.length() - 2)) && is_cytosine(s[i]) && !is_guanine(s[i + 1]) &&
         !is_guanine(s[i + 2]);
}

static inline bool
is_ddg(const std::string &s, std::size_t i) {
  return (i < (s.length() - 2)) && !is_cytosine(s[i]) &&
         !is_cytosine(s[i + 1]) && is_guanine(s[i + 2]);
}

static inline bool
is_c_at_g(const std::string &s, std::size_t i) {
  return (i < (s.length() - 2)) && is_cytosine(s[i]) &&
         !is_cytosine(s[i + 1]) && !is_guanine(s[i + 1]) &&
         is_guanine(s[i + 2]);
}

/* The "tag" returned by this function should be exclusive, so that
 * the order of checking conditions doesn't matter. There is also a
 * bit of a hack in that the unsigned "pos" could wrap, but this still
 * works as long as the chromosome size is not the maximum size of a
 * std::size_t.
 */
static inline std::uint32_t
get_tag_from_genome_c(const std::string &s, const std::size_t pos) {
  if (is_cpg(s, pos))
    return 0;
  else if (is_chh(s, pos))
    return 1;
  else if (is_c_at_g(s, pos))
    return 2;
  return 3;
}

static inline std::uint32_t
get_tag_from_genome_g(const std::string &s, const std::size_t pos) {
  if (is_cpg(s, pos - 1))
    return 0;
  else if (is_ddg(s, pos - 2))
    return 1;
  else if (is_c_at_g(s, pos - 2))
    return 2;
  return 3;
}

static const char *tag_values[] = {
  "CpG",  // 0
  "CHH",  // 1
  "CXG",  // 2
  "CCG",  // 3
  "N"     // 4
};

static void
write_missing_sites(const std::string &name, const std::string &chrom,
                    const std::uint64_t start_pos, const std::uint64_t end_pos,
                    bamxx::bgzf_file &out) {
  const std::string name_tab = name + "\t";
  quick_buf buf;
  for (auto pos = start_pos; pos < end_pos; ++pos) {
    const char base = chrom[pos];
    if (is_cytosine(base) || is_guanine(base)) {
      const bool is_c = is_cytosine(base);
      const std::uint32_t the_tag = is_c ? get_tag_from_genome_c(chrom, pos)
                                         : get_tag_from_genome_g(chrom, pos);
      buf.clear();
      buf << name_tab << pos << (is_c ? "\t+\t" : "\t-\t")
          << tag_values[the_tag] << "\t0\t0\n";
      if (!out.write(buf.c_str(), buf.tellp()))
        throw std::runtime_error("error writing output");
    }
  }
}

static void
write_current_site(const MSite &site, bamxx::bgzf_file &out) {
  quick_buf buf;  // keep underlying buffer space?
  buf << site << '\n';
  if (!out.write(buf.c_str(), buf.tellp()))
    throw std::runtime_error("error writing site: " + site.tostring());
}

typedef std::vector<std::string>::const_iterator chrom_itr_t;

static chrom_itr_t
get_chrom(const std::unordered_map<std::string, chrom_itr_t> &chrom_lookup,
          const std::string &chrom_name) {
  const auto chrom_idx = chrom_lookup.find(chrom_name);
  if (chrom_idx == cend(chrom_lookup))
    throw std::runtime_error("chromosome not found: " + chrom_name);
  return chrom_idx->second;
}

static std::int32_t
get_chrom_idx(const std::unordered_map<std::string, std::int32_t> &name_to_idx,
              const std::string &chrom_name) {
  const auto chrom_idx = name_to_idx.find(chrom_name);
  if (chrom_idx == cend(name_to_idx))
    throw std::runtime_error("chromosome not found: " + chrom_name);
  return chrom_idx->second;
}

static void
process_sites(const bool verbose, const bool add_missing_chroms,
              const bool compress_output, const std::size_t n_threads,
              const std::string &infile, const std::string &outfile,
              const std::string &chroms_file) {

  // first get the chromosome names and sequences from the FASTA file
  std::vector<std::string> chroms, names;
  read_fasta_file_short_names(chroms_file, names, chroms);
  for (auto &i : chroms)
    transform(cbegin(i), cend(i), begin(i),
              [](const char c) { return std::toupper(c); });
  if (verbose)
    std::cerr << "[n chroms in reference: " << chroms.size() << "]\n";

  std::unordered_map<std::string, chrom_itr_t> chrom_lookup;
  std::unordered_map<std::string, std::int32_t> name_to_idx;
  std::vector<std::uint64_t> chrom_sizes(size(chroms), 0);
  for (std::size_t i = 0; i < size(chroms); ++i) {
    chrom_lookup[names[i]] = cbegin(chroms) + i;
    name_to_idx[names[i]] = i;
    chrom_sizes[i] = size(chroms[i]);
  }

  if (add_missing_chroms)
    verify_chrom_orders(verbose, n_threads, infile, name_to_idx);

  bamxx::bam_tpool tp(n_threads);

  bamxx::bgzf_file in(infile, "r");
  if (!in)
    throw std::runtime_error("failed to open input file");

  const std::string output_mode = compress_output ? "w" : "wu";
  bamxx::bgzf_file out(outfile, output_mode);
  if (!out)
    throw std::runtime_error("error opening output file: " + outfile);

  // set the threads for the input file decompression
  if (n_threads > 1) {
    tp.set_io(in);
    tp.set_io(out);
  }

  MSite site;
  std::string chrom_name;
  std::int32_t prev_chrom_idx = -1;
  std::uint64_t pos = std::numeric_limits<std::uint64_t>::max();

  // ADS: this is probably a poor strategy since we already would know
  // the index of the chrom sequence in the vector.
  chrom_itr_t chrom_itr;
  std::string line;

  while (getline(in, line)) {
    if (is_counts_header_line(line)) {
      write_counts_header_line(line, out);
      continue;
    }
    site.initialize(line.data(), line.data() + size(line));
    if (site.chrom != chrom_name) {

      if (pos != std::numeric_limits<std::uint64_t>::max())
        write_missing_sites(chrom_name, *chrom_itr, pos, size(*chrom_itr), out);

      const std::int32_t chrom_idx = get_chrom_idx(name_to_idx, site.chrom);

      if (add_missing_chroms)
        for (auto i = prev_chrom_idx + 1; i < chrom_idx; ++i) {
          if (verbose)
            std::cerr << "processing: " << names[i] << " (missing)\n";
          write_missing_sites(names[i], chroms[i], 0u, size(chroms[i]), out);
        }

      chrom_name = site.chrom;
      chrom_itr = get_chrom(chrom_lookup, chrom_name);
      pos = 0;
      prev_chrom_idx = chrom_idx;
      if (verbose)
        std::cerr << "processing: " << chrom_name << "\n";
    }
    if (pos < site.pos)
      write_missing_sites(chrom_name, *chrom_itr, pos, site.pos, out);
    write_current_site(site, out);
    pos = site.pos + 1;
  }
  write_missing_sites(chrom_name, *chrom_itr, pos, size(*chrom_itr), out);

  if (add_missing_chroms) {
    const std::int32_t chrom_idx = size(chroms);
    for (auto i = prev_chrom_idx + 1; i < chrom_idx; ++i) {
      if (verbose)
        std::cerr << "processing: " << names[i] << " (missing)\n";
      write_missing_sites(names[i], chroms[i], 0u, size(chroms[i]), out);
    }
  }
}

int
main_recovered(int argc, char *argv[]) {
  try {

    bool verbose = false;
    bool add_missing_chroms = false;
    bool compress_output = false;
    std::size_t n_threads = 1;

    std::string outfile;
    std::string chroms_file;
    const std::string description =
      "add sites that are missing as non-covered sites";

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), description,
                           "<methcounts-file>");
    opt_parse.add_opt("output", 'o', "output file (required)", true, outfile);
    opt_parse.add_opt("missing", 'm', "add missing chroms", false,
                      add_missing_chroms);
    opt_parse.add_opt("threads", 't', "number of threads", false, n_threads);
    opt_parse.add_opt("chrom", 'c', "reference genome file (FASTA format)",
                      true, chroms_file);
    opt_parse.add_opt("zip", 'z', "output gzip format", false, compress_output);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, verbose);
    std::vector<std::string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      std::cerr << opt_parse.help_message() << "\n"
                << opt_parse.about_message() << "\n";
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      std::cerr << opt_parse.about_message() << "\n";
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      std::cerr << opt_parse.option_missing_message() << "\n";
      return EXIT_SUCCESS;
    }
    if (leftover_args.size() != 1) {
      std::cerr << opt_parse.help_message() << "\n";
      return EXIT_SUCCESS;
    }
    const std::string filename(leftover_args.front());
    /****************** END COMMAND LINE OPTIONS *****************/

    process_sites(verbose, add_missing_chroms, compress_output, n_threads,
                  filename, outfile, chroms_file);
  }
  catch (const std::exception &e) {
    std::cerr << e.what() << "\n";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
