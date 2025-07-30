/* unxcounts: reverse the process of xcounts and generate the counts
 * file, including sites not covered.
 *
 * Copyright (C) 2023 Andrew D. Smith
 *
 * Authors: Andrew D. Smith
 *
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 */

#include "MSite.hpp"
#include "bsutils.hpp"
#include "counts_header.hpp"

#include <bamxx.hpp>

// from smithlab_cpp
#include "OptionParser.hpp"
#include "smithlab_os.hpp"
#include "smithlab_utils.hpp"

#include <charconv>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

static void
read_fasta_file_short_names_uppercase(const std::string &chroms_file,
                                      std::vector<std::string> &names,
                                      std::vector<std::string> &chroms) {
  chroms.clear();
  names.clear();
  read_fasta_file_short_names(chroms_file, names, chroms);
  for (auto &i : chroms)
    transform(std::cbegin(i), std::cend(i), begin(i),
              [](const char c) { return std::toupper(c); });
}

static void
verify_chrom_orders(
  const bool verbose, const std::uint32_t n_threads,
  const std::string &filename,
  const std::unordered_map<std::string, std::int32_t> &chroms_order) {
  bamxx::bam_tpool tp(n_threads);

  bamxx::bgzf_file in(filename, "r");
  if (!in)
    throw std::runtime_error("bad file: " + filename);

  // set the threads for the input file decompression
  if (n_threads > 1 && in.is_bgzf())
    tp.set_io(in);

  std::unordered_set<std::int32_t> chroms_seen;
  std::int32_t prev_id = -1;

  kstring_t line = KS_INITIALIZE;
  const int ret = ks_resize(&line, 1024);
  if (ret)
    throw std::runtime_error("failed to acquire buffer");

  while (bamxx::getline(in, line)) {
    if (std::isdigit(line.s[0]))
      continue;
    if (is_counts_header_line(line.s))
      continue;

    std::string chrom{line.s};
    if (verbose)
      std::cerr << "verifying: " << chrom << "\n";

    const auto idx_itr = chroms_order.find(chrom);
    if (idx_itr == std::cend(chroms_order))
      throw std::runtime_error("chrom not found genome file: " + chrom);
    const auto idx = idx_itr->second;

    if (chroms_seen.find(idx) != end(chroms_seen))
      throw std::runtime_error("chroms out of order in: " + filename);
    chroms_seen.insert(idx);

    if (idx < prev_id)
      throw std::runtime_error("inconsistent chromosome order at: " + chrom);

    prev_id = idx;
  }
  if (verbose)
    std::cerr << "chrom orders are consistent" << "\n";
}

static const char *tag_values[] = {
  "CpG",  // 0
  "CHH",  // 1
  "CXG",  // 2
  "CCG",  // 3
  "N"     // 4
};

static const int tag_sizes[] = {3, 3, 3, 3, 1};

// ADS: the values below allow for things like CHH where the is a N in
// the triplet; I'm allowing that for consistency with the weird logic
// from earlier versions.
const std::uint32_t context_codes[] = {
  /*CAA CHH*/ 1,
  /*CAC CHH*/ 1,
  /*CAG CXG*/ 2,
  /*CAT CHH*/ 1,
  /*CAN ---*/ 1,  // 4,
  /*CCA CHH*/ 1,
  /*CCC CHH*/ 1,
  /*CCG CCG*/ 3,
  /*CCT CHH*/ 1,
  /*CCN ---*/ 1,  // 4,
  /*CGA CpG*/ 0,
  /*CGC CpG*/ 0,
  /*CGG CpG*/ 0,
  /*CGT CpG*/ 0,
  /*CGN CpG*/ 0,
  /*CTA CHH*/ 1,
  /*CTC CHH*/ 1,
  /*CTG CXG*/ 2,
  /*CTT CHH*/ 1,
  /*CTN ---*/ 1,  // 4,
  /*CNA ---*/ 1,  // 4,
  /*CNC ---*/ 1,  // 4,
  /*CNG ---*/ 2,
  /*CNT ---*/ 1,  // 4,
  /*CNN ---*/ 1   // 4
};

static inline std::uint32_t
get_tag_from_genome_c(const std::string &s, const size_t pos) {
  const auto val = base2int(s[pos + 1]) * 5 + base2int(s[pos + 2]);
  return context_codes[val];
}

static inline std::uint32_t
get_tag_from_genome_g(const std::string &s, const size_t pos) {
  const auto val =
    base2int(complement(s[pos - 1])) * 5 + base2int(complement(s[pos - 2]));
  return context_codes[val];
}

static bool
write_missing(const std::uint32_t name_size, const std::string &chrom,
              const std::uint64_t start_pos, const std::uint64_t end_pos,
              std::vector<char> &buf, bamxx::bgzf_file &out) {
  static constexpr auto zeros = "\t0\t0\n";
  static constexpr auto pos_strand = "\t+\t";
  static constexpr auto neg_strand = "\t-\t";
  const auto buf_end = buf.data() + size(buf);
  // chrom name is already in the buffer so move past it
  auto cursor = buf.data() + name_size + 1;
  for (auto pos = start_pos; pos < end_pos; ++pos) {
    const char base = chrom[pos];
    if (is_cytosine(base) || is_guanine(base)) {
      const bool is_c = is_cytosine(base);
      const std::uint32_t the_tag = is_c ? get_tag_from_genome_c(chrom, pos)
                                         : get_tag_from_genome_g(chrom, pos);
#pragma GCC diagnostic push
#pragma GCC diagnostic error "-Wstringop-overflow=0"
      auto [ptr, ec] = std::to_chars(cursor, buf_end, pos);
      ptr = std::copy_n(is_c ? pos_strand : neg_strand, 3, ptr);
      ptr = std::copy_n(tag_values[the_tag], tag_sizes[the_tag], ptr);
      ptr = std::copy_n(zeros, 5, ptr);
      const auto sz = std::distance(buf.data(), ptr);
#pragma GCC diagnostic push

      if (bgzf_write(out.f, buf.data(), sz) != sz)
        return false;
    }
  }
  return true;
}

static bool
write_missing_cpg(const std::uint32_t &name_size, const std::string &chrom,
                  const std::uint64_t start_pos, const std::uint64_t end_pos,
                  std::vector<char> &buf, bamxx::bgzf_file &out) {
  static constexpr auto zeros = "\t0\t0\n";
  static constexpr auto pos_strand = "\t+\t";
  const auto buf_end = buf.data() + size(buf);
  // chrom name is already in the buffer so move past it
  auto cursor = buf.data() + name_size + 1;
  for (auto pos = start_pos; pos < end_pos - 1; ++pos) {
    // When this function is called, the "end_pos" is either the chrom
    // size or the position of a base known to be a C. So we never
    // have to allow pos+1 to equal end_pos.
    if (is_cytosine(chrom[pos]) && is_guanine(chrom[pos + 1])) {
#pragma GCC diagnostic push
#pragma GCC diagnostic error "-Wstringop-overflow=0"
      auto [ptr, ec] = std::to_chars(cursor, buf_end, pos);
      ptr = std::copy_n(pos_strand, 3, ptr);
      ptr = std::copy_n("CpG", 3, ptr);
      ptr = std::copy_n(zeros, 5, ptr);
      const auto sz = std::distance(buf.data(), ptr);
#pragma GCC diagnostic push
      if (bgzf_write(out.f, buf.data(), sz) != sz)
        return false;
    }
  }
  return true;
}

static bool
write_site(const std::uint32_t name_size, const std::string &chrom,
           const std::uint32_t pos, const std::uint32_t n_meth,
           const std::uint32_t n_unmeth, std::vector<char> &buf,
           bamxx::bgzf_file &out) {
  static constexpr auto pos_strand = "\t+\t";
  static constexpr auto neg_strand = "\t-\t";
  static constexpr auto fmt = std::chars_format::general;

  const auto buf_end = buf.data() + size(buf);
  const char base = chrom[pos];
  assert(is_cytosine(base) || is_guanine(base));
  const bool is_c = is_cytosine(base);
  const std::uint32_t the_tag = is_c ? get_tag_from_genome_c(chrom, pos)
                                     : get_tag_from_genome_g(chrom, pos);
  const auto n_reads = n_meth + n_unmeth;
  const auto meth = static_cast<double>(n_meth) / std::max(n_reads, 1u);

#pragma GCC diagnostic push
#pragma GCC diagnostic error "-Wstringop-overflow=0"
  // chrom name is already in the buffer so move past it
  auto cursor = buf.data() + name_size + 1;
  {
    auto [ptr, ec] = std::to_chars(cursor, buf_end, pos);
    cursor = ptr;
  }
  cursor = std::copy_n(is_c ? pos_strand : neg_strand, 3, cursor);
  cursor = std::copy_n(tag_values[the_tag], tag_sizes[the_tag], cursor);
  *cursor++ = '\t';
  {
    // use default precision, 6, same as std::cout default
    auto [ptr, ec] = std::to_chars(cursor, buf_end, meth, fmt, 6);
    cursor = ptr;
  }
  *cursor++ = '\t';
  {
    auto [ptr, ec] = std::to_chars(cursor, buf_end, n_reads);
    cursor = ptr;
  }
  *cursor++ = '\n';
#pragma GCC diagnostic push

  const auto sz = std::distance(buf.data(), cursor);
  return bgzf_write(out.f, buf.data(), sz) == sz;
}

typedef std::vector<std::string>::const_iterator chrom_itr_t;

static chrom_itr_t
get_chrom(const std::unordered_map<std::string, chrom_itr_t> &chrom_lookup,
          const std::string &chrom_name) {
  const auto chr_id = chrom_lookup.find(chrom_name);
  if (chr_id == std::cend(chrom_lookup))
    throw std::runtime_error("chromosome not found: " + chrom_name);
  return chr_id->second;
}

static std::int32_t
get_chrom_id(const std::unordered_map<std::string, std::int32_t> &name_to_id,
             const std::string &chrom_name) {
  const auto chr_id = name_to_id.find(chrom_name);
  if (chr_id == std::cend(name_to_id))
    throw std::runtime_error("chromosome not found: " + chrom_name);
  return chr_id->second;
}

static bool
verify_chrom(const std::string &header_line,
             const std::unordered_map<std::string, std::int32_t> &name_to_id,
             const std::vector<std::uint64_t> &chrom_sizes) {
  if (is_counts_header_version_line(header_line))
    return true;
  std::istringstream iss(header_line.substr(1));
  std::string name;
  std::uint64_t chrom_size = 0;
  if (!(iss >> name >> chrom_size))
    return false;

  const auto idx = name_to_id.find(name);
  if (idx == std::cend(name_to_id))
    return false;

  return chrom_size == chrom_sizes[idx->second];
}

static void
get_lookups(const std::vector<std::string> &names,
            const std::vector<std::string> &chroms,
            std::unordered_map<std::string, chrom_itr_t> &chrom_lookup,
            std::unordered_map<std::string, std::int32_t> &name_to_id,
            std::vector<std::uint64_t> &chrom_sizes) {
  chrom_lookup.clear();
  name_to_id.clear();
  chrom_sizes = std::vector<std::uint64_t>(size(chroms), 0);
  for (size_t i = 0; i < size(chroms); ++i) {
    chrom_lookup[names[i]] = std::cbegin(chroms) + i;
    name_to_id[names[i]] = i;
    chrom_sizes[i] = size(chroms[i]);
  }
}

static void
process_header_line(
  const std::unordered_map<std::string, std::int32_t> &name_to_id,
  const std::vector<std::uint64_t> &chrom_sizes, const kstring_t &line,
  bamxx::bgzf_file &out) {
  std::string hdr_line{line.s};
  if (size(hdr_line) > 1 && !verify_chrom(hdr_line, name_to_id, chrom_sizes))
    throw std::runtime_error{"failed to verify header for: " + hdr_line};
  if (!write_counts_header_line(hdr_line, out))
    throw std::runtime_error{"failed to write header line: " + hdr_line};
}

// write all sites for chroms in the given range
static void
write_all_sites(const bool verbose, const std::uint32_t prev_chr_id,
                const std::uint32_t chr_id,
                const std::vector<std::string> &names,
                const std::vector<std::string> &chroms, std::vector<char> &buf,
                bamxx::bgzf_file &out) {
  for (auto i = prev_chr_id + 1; i < chr_id; ++i) {
    if (verbose)
      std::cerr << "processing: " << names[i] << " (missing)" << "\n";
    auto res =
      std::copy(std::cbegin(names[i]), std::cend(names[i]), buf.data());
    *res = '\t';
    write_missing(size(names[i]), chroms[i], 0u, size(chroms[i]), buf, out);
  }
}

static void
process_sites(const bool verbose, const bool add_missing_chroms,
              const bool require_covered, const bool compress_output,
              const size_t n_threads, const std::string &infile,
              const std::string &outfile, const std::string &chroms_file) {
  // first get the chromosome names and sequences from the FASTA file
  std::vector<std::string> chroms, names;
  read_fasta_file_short_names_uppercase(chroms_file, names, chroms);
  if (verbose)
    std::cerr << "[n chroms in reference: " << chroms.size() << "]" << "\n";

  std::unordered_map<std::string, chrom_itr_t> chrom_lookup;
  std::unordered_map<std::string, std::int32_t> name_to_id;
  std::vector<std::uint64_t> chrom_sizes(size(chroms), 0);
  get_lookups(names, chroms, chrom_lookup, name_to_id, chrom_sizes);

  if (add_missing_chroms)
    verify_chrom_orders(verbose, n_threads, infile, name_to_id);

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
    if (in.is_bgzf())
      tp.set_io(in);
    tp.set_io(out);
  }

  static constexpr std::uint32_t output_buffer_size = 1024;
  std::vector<char> buf(output_buffer_size, '\0');

  kstring_t line = KS_INITIALIZE;
  const int ret = ks_resize(&line, output_buffer_size);
  if (ret)
    throw std::runtime_error("failed to acquire buffer");

  std::string chrom_name;
  std::uint32_t nm_sz{};
  std::int32_t prev_chr_id = -1;
  std::uint64_t pos = std::numeric_limits<std::uint64_t>::max();

  // ADS: this is probably a poor strategy since we already would know
  // the index of the chrom sequence in the vector.
  chrom_itr_t ch_itr;

  while (getline(in, line)) {
    if (is_counts_header_line(line.s)) {
      process_header_line(name_to_id, chrom_sizes, line, out);
      continue;  // ADS: early loop exit
    }

    if (!std::isdigit(line.s[0])) {  // check if we have a chrom line

      if (!require_covered && pos != std::numeric_limits<std::uint64_t>::max())
        write_missing(nm_sz, *ch_itr, pos + 1, size(*ch_itr), buf, out);

      chrom_name = std::string{line.s};
      nm_sz = size(chrom_name);
      const std::int32_t chr_id = get_chrom_id(name_to_id, chrom_name);

      if (add_missing_chroms)
        write_all_sites(verbose, prev_chr_id, chr_id, names, chroms, buf, out);

      ch_itr = get_chrom(chrom_lookup, chrom_name);
      pos = 0;
      prev_chr_id = chr_id;
      if (verbose)
        std::cerr << "processing: " << chrom_name << "\n";

      auto res =
        std::copy(std::cbegin(chrom_name), std::cend(chrom_name), buf.data());
      *res = '\t';
    }
    else {
      std::uint32_t pos_step = 0, n_meth = 0, n_unmeth = 0;
      const auto end_line = line.s + line.l;
      auto res = std::from_chars(line.s, end_line, pos_step);
      res = std::from_chars(res.ptr + 1, end_line, n_meth);
      res = std::from_chars(res.ptr + 1, end_line, n_unmeth);

      const auto curr_pos = pos + pos_step;
      if (!require_covered && pos + 1 < curr_pos)
        write_missing(nm_sz, *ch_itr, pos + 1, curr_pos, buf, out);

      write_site(nm_sz, *ch_itr, curr_pos, n_meth, n_unmeth, buf, out);
      pos = curr_pos;
    }
  }
  if (!require_covered)
    write_missing(nm_sz, *ch_itr, pos + 1, size(*ch_itr), buf, out);
  if (add_missing_chroms)
    write_all_sites(verbose, prev_chr_id, size(chroms), names, chroms, buf,
                    out);
}

// write all cpg sites for chroms in the given range
static void
write_all_cpgs(const bool verbose, const std::uint32_t prev_chr_id,
               const std::uint32_t chr_id,
               const std::vector<std::string> &names,
               const std::vector<std::string> &chroms, std::vector<char> &buf,
               bamxx::bgzf_file &out) {
  for (auto i = prev_chr_id + 1; i < chr_id; ++i) {
    if (verbose)
      std::cerr << "processing: " << names[i] << " (missing)" << "\n";
    auto res =
      std::copy(std::cbegin(names[i]), std::cend(names[i]), buf.data());
    *res = '\t';
    write_missing_cpg(size(names[i]), chroms[i], 0u, size(chroms[i]), buf, out);
  }
}

static void
process_cpg_sites(const bool verbose, const bool add_missing_chroms,
                  const bool require_covered, const bool compress_output,
                  const size_t n_threads, const std::string &infile,
                  const std::string &outfile, const std::string &chroms_file) {
  // first get the chromosome names and sequences from the FASTA file
  std::vector<std::string> chroms, names;
  read_fasta_file_short_names_uppercase(chroms_file, names, chroms);
  if (verbose)
    std::cerr << "[n chroms in reference: " << chroms.size() << "]" << "\n";

  std::unordered_map<std::string, chrom_itr_t> chrom_lookup;
  std::unordered_map<std::string, std::int32_t> name_to_id;
  std::vector<std::uint64_t> chrom_sizes(size(chroms), 0);
  get_lookups(names, chroms, chrom_lookup, name_to_id, chrom_sizes);

  if (add_missing_chroms)
    verify_chrom_orders(verbose, n_threads, infile, name_to_id);

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
    if (in.is_bgzf())
      tp.set_io(in);
    tp.set_io(out);
  }

  static constexpr std::uint32_t output_buffer_size = 1024;
  std::vector<char> buf(output_buffer_size, '\0');

  kstring_t line = KS_INITIALIZE;
  const int ret = ks_resize(&line, output_buffer_size);
  if (ret)
    throw std::runtime_error("failed to acquire buffer");

  std::string chrom_name;
  std::uint32_t nm_sz{};
  std::int32_t prev_chr_id = -1;
  std::uint64_t pos = std::numeric_limits<std::uint64_t>::max();

  // ADS: this is probably a poor strategy since we already would know
  // the index of the chrom sequence in the vector.
  chrom_itr_t ch_itr;

  while (getline(in, line)) {
    if (is_counts_header_line(line.s)) {
      process_header_line(name_to_id, chrom_sizes, line, out);
      continue;  // ADS: early loop exit
    }

    if (!std::isdigit(line.s[0])) {  // check if we have a chrom line

      if (!require_covered && pos != std::numeric_limits<std::uint64_t>::max())
        write_missing_cpg(nm_sz, *ch_itr, pos + 1, size(*ch_itr), buf, out);

      chrom_name = std::string{line.s};
      nm_sz = size(chrom_name);
      const std::int32_t chr_id = get_chrom_id(name_to_id, chrom_name);

      if (add_missing_chroms)
        write_all_cpgs(verbose, prev_chr_id, chr_id, names, chroms, buf, out);

      ch_itr = get_chrom(chrom_lookup, chrom_name);
      pos = 0;
      prev_chr_id = chr_id;
      if (verbose)
        std::cerr << "processing: " << chrom_name << "\n";

      auto res =
        std::copy(std::cbegin(chrom_name), std::cend(chrom_name), buf.data());
      *res = '\t';
    }
    else {
      std::uint32_t pos_step = 0, n_meth = 0, n_unmeth = 0;
      const auto end_line = line.s + line.l;
      auto res = std::from_chars(line.s, end_line, pos_step);
      res = std::from_chars(res.ptr + 1, end_line, n_meth);
      res = std::from_chars(res.ptr + 1, end_line, n_unmeth);

      const auto curr_pos = pos + pos_step;
      if (!require_covered && pos + 1 < curr_pos)
        write_missing_cpg(nm_sz, *ch_itr, pos + 1, curr_pos, buf, out);

      write_site(nm_sz, *ch_itr, curr_pos, n_meth, n_unmeth, buf, out);
      pos = curr_pos;
    }
  }
  if (!require_covered)
    write_missing_cpg(nm_sz, *ch_itr, pos + 1, size(*ch_itr), buf, out);
  if (add_missing_chroms)
    write_all_cpgs(verbose, prev_chr_id, size(chroms), names, chroms, buf, out);
}

int
main_unxcounts(int argc, char *argv[]) {
  try {
    bool verbose = false;
    bool add_missing_chroms = false;
    bool require_covered = false;
    bool compress_output = false;
    bool assume_cpg_only = false;
    size_t n_threads = 1;

    std::string outfile;
    std::string chroms_file;
    const std::string description =
      "convert compressed counts format back to full counts";

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), description, "<xcounts-file>");
    opt_parse.add_opt("output", 'o', "output file (required)", true, outfile);
    opt_parse.add_opt("missing", 'm', "add missing chroms", false,
                      add_missing_chroms);
    opt_parse.add_opt("cpg", '\0', "assume symmetric CpG input and output",
                      false, assume_cpg_only);
    opt_parse.add_opt("reads", 'r', "report only sites with reads", false,
                      require_covered);
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
    if (require_covered && add_missing_chroms) {
      std::cerr << "options mutually exclusive: reads and missing" << "\n";
      return EXIT_FAILURE;
    }
    const std::string filename(leftover_args.front());
    /****************** END COMMAND LINE OPTIONS *****************/

    if (require_covered && add_missing_chroms) {
      std::cerr << "options mutually exclusive: reads and missing" << "\n";
      return EXIT_FAILURE;
    }

    if (assume_cpg_only)
      process_cpg_sites(verbose, add_missing_chroms, require_covered,
                        compress_output, n_threads, filename, outfile,
                        chroms_file);
    else
      process_sites(verbose, add_missing_chroms, require_covered,
                    compress_output, n_threads, filename, outfile, chroms_file);
  }
  catch (const std::exception &e) {
    std::cerr << e.what() << "\n";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
