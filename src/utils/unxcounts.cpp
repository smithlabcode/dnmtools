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

#include <bamxx.hpp>

#include <charconv>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// from smithlab_cpp
#include "MSite.hpp"
#include "OptionParser.hpp"
#include "bsutils.hpp"
#include "counts_header.hpp"
#include "dnmt_error.hpp"
#include "smithlab_os.hpp"
#include "smithlab_utils.hpp"

using std::cbegin;
using std::cend;
using std::cerr;
using std::copy;
using std::copy_n;
using std::cout;
using std::endl;
using std::from_chars;
using std::numeric_limits;
using std::pair;
using std::runtime_error;
using std::string;
using std::to_chars;
using std::to_string;
using std::unordered_map;
using std::unordered_set;
using std::vector;

using bamxx::bgzf_file;

template<typename T> using num_lim = std::numeric_limits<T>;

static void
read_fasta_file_short_names_uppercase(const string &chroms_file,
                                      vector<string> &names,
                                      vector<string> &chroms) {
  chroms.clear();
  names.clear();
  read_fasta_file_short_names(chroms_file, names, chroms);
  for (auto &i : chroms)
    transform(cbegin(i), cend(i), begin(i),
              [](const char c) { return std::toupper(c); });
}


static void
verify_chrom_orders(const bool verbose, const uint32_t n_threads,
                    const string &filename,
                    const unordered_map<string, int32_t> &chroms_order) {
  bamxx::bam_tpool tp(n_threads);

  bgzf_file in(filename, "r");
  if (!in) throw runtime_error("bad file: " + filename);

  // set the threads for the input file decompression
  if (n_threads > 1 && in.is_bgzf()) tp.set_io(in);

  unordered_set<int32_t> chroms_seen;
  int32_t prev_id = -1;

  kstring_t line{0, 0, nullptr};
  const int ret = ks_resize(&line, 1024);
  if (ret) throw runtime_error("failed to acquire buffer");

  while (bamxx::getline(in, line)) {
    if (std::isdigit(line.s[0])) continue;
    if (is_counts_header_line(line.s)) continue;

    string chrom{line.s};
    if (verbose) cerr << "verifying: " << chrom << endl;

    const auto idx_itr = chroms_order.find(chrom);
    if (idx_itr == cend(chroms_order))
      throw runtime_error("chrom not found genome file: " + chrom);
    const auto idx = idx_itr->second;

    if (chroms_seen.find(idx) != end(chroms_seen))
      throw runtime_error("chroms out of order in: " + filename);
    chroms_seen.insert(idx);

    if (idx < prev_id)
      throw runtime_error("inconsistent chromosome order at: " + chrom);

    prev_id = idx;
  }
  if (verbose) cerr << "chrom orders are consistent" << endl;
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
const uint32_t context_codes[] = {
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

static inline uint32_t
get_tag_from_genome_c(const string &s, const size_t pos) {
  const auto val = base2int(s[pos + 1]) * 5 + base2int(s[pos + 2]);
  return context_codes[val];
}

static inline uint32_t
get_tag_from_genome_g(const string &s, const size_t pos) {
  const auto val =
    base2int(complement(s[pos - 1])) * 5 + base2int(complement(s[pos - 2]));
  return context_codes[val];
}

static bool
write_missing(const uint32_t name_size, const string &chrom,
              const uint64_t start_pos, const uint64_t end_pos,
              vector<char> &buf, bgzf_file &out) {
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
      const uint32_t the_tag = is_c ? get_tag_from_genome_c(chrom, pos)
                                    : get_tag_from_genome_g(chrom, pos);
#pragma GCC diagnostic push
#pragma GCC diagnostic error "-Wstringop-overflow=0"
      auto [ptr, ec] = to_chars(cursor, buf_end, pos);
      ptr = copy_n(is_c ? pos_strand : neg_strand, 3, ptr);
      ptr = copy_n(tag_values[the_tag], tag_sizes[the_tag], ptr);
      ptr = copy_n(zeros, 5, ptr);
      const auto sz = std::distance(buf.data(), ptr);
#pragma GCC diagnostic push

      if (bgzf_write(out.f, buf.data(), sz) != sz) return false;
    }
  }
  return true;
}

static bool
write_missing_cpg(const uint32_t &name_size, const string &chrom,
                  const uint64_t start_pos, const uint64_t end_pos,
                  vector<char> &buf, bgzf_file &out) {
  static constexpr auto zeros = "\t0\t0\n";
  static constexpr auto pos_strand = "\t+\t";
  const auto buf_end = buf.data() + size(buf);
  // chrom name is already in the buffer so move past it
  auto cursor = buf.data() + name_size + 1;
  for (auto pos = start_pos; pos < end_pos - 1; ++pos) {
    // When this function is called, the "end_pos" is either the chrom
    // size or the position of a base known to be a C. So we never
    // have to allow pos+1 to equal end_pos.
    if (is_cytosine(chrom[pos]) && is_guanine(chrom[pos+1])) {
#pragma GCC diagnostic push
#pragma GCC diagnostic error "-Wstringop-overflow=0"
      auto [ptr, ec] = to_chars(cursor, buf_end, pos);
      ptr = copy_n(pos_strand, 3, ptr);
      ptr = copy_n("CpG", 3, ptr);
      ptr = copy_n(zeros, 5, ptr);
      const auto sz = std::distance(buf.data(), ptr);
#pragma GCC diagnostic push
      if (bgzf_write(out.f, buf.data(), sz) != sz) return false;
    }
  }
  return true;
}

static bool
write_site(const uint32_t name_size, const string &chrom, const uint32_t pos,
           const uint32_t n_meth, const uint32_t n_unmeth, vector<char> &buf,
           bgzf_file &out) {
  static constexpr auto pos_strand = "\t+\t";
  static constexpr auto neg_strand = "\t-\t";
  static constexpr auto fmt = std::chars_format::general;

  const auto buf_end = buf.data() + size(buf);
  const char base = chrom[pos];
  assert(is_cytosine(base) || is_guanine(base));
  const bool is_c = is_cytosine(base);
  const uint32_t the_tag = is_c ? get_tag_from_genome_c(chrom, pos)
                                : get_tag_from_genome_g(chrom, pos);
  const auto n_reads = n_meth + n_unmeth;
  const auto meth = static_cast<double>(n_meth) / std::max(n_reads, 1u);

#pragma GCC diagnostic push
#pragma GCC diagnostic error "-Wstringop-overflow=0"
  // chrom name is already in the buffer so move past it
  auto cursor = buf.data() + name_size + 1;
  {
    auto [ptr, ec] = to_chars(cursor, buf_end, pos);
    cursor = ptr;
  }
  cursor = copy_n(is_c ? pos_strand : neg_strand, 3, cursor);
  cursor = copy_n(tag_values[the_tag], tag_sizes[the_tag], cursor);
  *cursor++ = '\t';
  {
    // use default precision, 6, same as cout default
    auto [ptr, ec] = to_chars(cursor, buf_end, meth, fmt, 6);
    cursor = ptr;
  }
  *cursor++ = '\t';
  {
    auto [ptr, ec] = to_chars(cursor, buf_end, n_reads);
    cursor = ptr;
  }
  *cursor++ = '\n';
#pragma GCC diagnostic push

  const auto sz = std::distance(buf.data(), cursor);
  return bgzf_write(out.f, buf.data(), sz) == sz;
}

typedef vector<string>::const_iterator chrom_itr_t;

static chrom_itr_t
get_chrom(const unordered_map<string, chrom_itr_t> &chrom_lookup,
          const string &chrom_name) {
  const auto chr_id = chrom_lookup.find(chrom_name);
  if (chr_id == cend(chrom_lookup))
    throw dnmt_error("chromosome not found: " + chrom_name);
  return chr_id->second;
}

static int32_t
get_chrom_id(const unordered_map<string, int32_t> &name_to_id,
              const string &chrom_name) {
  const auto chr_id = name_to_id.find(chrom_name);
  if (chr_id == cend(name_to_id))
    throw dnmt_error("chromosome not found: " + chrom_name);
  return chr_id->second;
}

static bool
verify_chrom(const string &header_line,
             const unordered_map<string, int32_t> &name_to_id,
             const vector<uint64_t> &chrom_sizes) {
  if (is_counts_header_version_line(header_line)) return true;
  std::istringstream iss(header_line.substr(1));
  string name;
  uint64_t chrom_size = 0;
  if (!(iss >> name >> chrom_size)) return false;

  const auto idx = name_to_id.find(name);
  if (idx == cend(name_to_id)) return false;

  return chrom_size == chrom_sizes[idx->second];
}

static void
get_lookups(const vector<string> &names, const vector<string> &chroms,
            unordered_map<string, chrom_itr_t> &chrom_lookup,
            unordered_map<string, int32_t> &name_to_id,
            vector<uint64_t> &chrom_sizes) {
  chrom_lookup.clear();
  name_to_id.clear();
  chrom_sizes = vector<uint64_t>(size(chroms), 0);
  for (size_t i = 0; i < size(chroms); ++i) {
    chrom_lookup[names[i]] = cbegin(chroms) + i;
    name_to_id[names[i]] = i;
    chrom_sizes[i] = size(chroms[i]);
  }
}

static void
process_header_line(const unordered_map<string, int32_t> &name_to_id,
                    const vector<uint64_t> &chrom_sizes, const kstring_t &line,
                    bgzf_file &out) {
  string hdr_line{line.s};
  if (size(hdr_line) > 1 && !verify_chrom(hdr_line, name_to_id, chrom_sizes))
    throw runtime_error{"failed to verify header for: " + hdr_line};
  if (!write_counts_header_line(hdr_line, out))
    throw runtime_error{"failed to write header line: " + hdr_line};
}


// write all sites for chroms in the given range
static void
write_all_sites(const bool verbose,
                const uint32_t prev_chr_id,
                const uint32_t chr_id,
                const vector<string> &names,
                const vector<string> &chroms,
                vector<char> &buf, bgzf_file &out) {
  for (auto i = prev_chr_id + 1; i < chr_id; ++i) {
    if (verbose)
      cerr << "processing: " << names[i] << " (missing)" << endl;
    auto res = copy(cbegin(names[i]), cend(names[i]), buf.data());
    *res = '\t';
    write_missing(size(names[i]), chroms[i], 0u, size(chroms[i]), buf, out);
  }
}

static void
process_sites(const bool verbose, const bool add_missing_chroms,
              const bool require_covered, const bool compress_output,
              const size_t n_threads, const string &infile,
              const string &outfile, const string &chroms_file) {
  // first get the chromosome names and sequences from the FASTA file
  vector<string> chroms, names;
  read_fasta_file_short_names_uppercase(chroms_file, names, chroms);
  if (verbose)
    cerr << "[n chroms in reference: " << chroms.size() << "]" << endl;

  unordered_map<string, chrom_itr_t> chrom_lookup;
  unordered_map<string, int32_t> name_to_id;
  vector<uint64_t> chrom_sizes(size(chroms), 0);
  get_lookups(names, chroms, chrom_lookup, name_to_id, chrom_sizes);

  if (add_missing_chroms)
    verify_chrom_orders(verbose, n_threads, infile, name_to_id);

  bamxx::bam_tpool tp(n_threads);

  bamxx::bgzf_file in(infile, "r");
  if (!in) throw dnmt_error("failed to open input file");

  const string output_mode = compress_output ? "w" : "wu";
  bgzf_file out(outfile, output_mode);
  if (!out) throw dnmt_error("error opening output file: " + outfile);

  // set the threads for the input file decompression
  if (n_threads > 1) {
    if (in.is_bgzf()) tp.set_io(in);
    tp.set_io(out);
  }

  static constexpr uint32_t output_buffer_size = 1024;
  vector<char> buf(output_buffer_size, '\0');

  kstring_t line{0, 0, nullptr};
  const int ret = ks_resize(&line, output_buffer_size);
  if (ret) throw runtime_error("failed to acquire buffer");

  string chrom_name;
  uint32_t nm_sz{};
  int32_t prev_chr_id = -1;
  uint64_t pos = num_lim<uint64_t>::max();

  // ADS: this is probably a poor strategy since we already would know
  // the index of the chrom sequence in the vector.
  chrom_itr_t ch_itr;

  while (getline(in, line)) {
    if (is_counts_header_line(line.s)) {
      process_header_line(name_to_id, chrom_sizes, line, out);
      continue;  // ADS: early loop exit
    }

    if (!std::isdigit(line.s[0])) {  // check if we have a chrom line

      if (!require_covered && pos != num_lim<uint64_t>::max())
        write_missing(nm_sz, *ch_itr, pos + 1, size(*ch_itr), buf, out);

      chrom_name = string{line.s};
      nm_sz = size(chrom_name);
      const int32_t chr_id = get_chrom_id(name_to_id, chrom_name);

      if (add_missing_chroms)
        write_all_sites(verbose, prev_chr_id, chr_id, names, chroms, buf, out);

      ch_itr = get_chrom(chrom_lookup, chrom_name);
      pos = 0;
      prev_chr_id = chr_id;
      if (verbose) cerr << "processing: " << chrom_name << endl;

      auto res = copy(cbegin(chrom_name), cend(chrom_name), buf.data());
      *res = '\t';
    }
    else {
      uint32_t pos_step = 0, n_meth = 0, n_unmeth = 0;
      const auto end_line = line.s + line.l;
      auto res = from_chars(line.s, end_line, pos_step);
      res = from_chars(res.ptr + 1, end_line, n_meth);
      res = from_chars(res.ptr + 1, end_line, n_unmeth);

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
write_all_cpgs(const bool verbose,
               const uint32_t prev_chr_id,
               const uint32_t chr_id,
               const vector<string> &names,
               const vector<string> &chroms,
               vector<char> &buf, bgzf_file &out) {
  for (auto i = prev_chr_id + 1; i < chr_id; ++i) {
    if (verbose)
      cerr << "processing: " << names[i] << " (missing)" << endl;
    auto res = copy(cbegin(names[i]), cend(names[i]), buf.data());
    *res = '\t';
    write_missing_cpg(size(names[i]), chroms[i], 0u, size(chroms[i]), buf, out);
  }
}

static void
process_cpg_sites(const bool verbose, const bool add_missing_chroms,
                  const bool require_covered, const bool compress_output,
                  const size_t n_threads, const string &infile,
                  const string &outfile, const string &chroms_file) {
  // first get the chromosome names and sequences from the FASTA file
  vector<string> chroms, names;
  read_fasta_file_short_names_uppercase(chroms_file, names, chroms);
  if (verbose)
    cerr << "[n chroms in reference: " << chroms.size() << "]" << endl;

  unordered_map<string, chrom_itr_t> chrom_lookup;
  unordered_map<string, int32_t> name_to_id;
  vector<uint64_t> chrom_sizes(size(chroms), 0);
  get_lookups(names, chroms, chrom_lookup, name_to_id, chrom_sizes);

  if (add_missing_chroms)
    verify_chrom_orders(verbose, n_threads, infile, name_to_id);

  bamxx::bam_tpool tp(n_threads);

  bamxx::bgzf_file in(infile, "r");
  if (!in) throw dnmt_error("failed to open input file");

  const string output_mode = compress_output ? "w" : "wu";
  bgzf_file out(outfile, output_mode);
  if (!out) throw dnmt_error("error opening output file: " + outfile);

  // set the threads for the input file decompression
  if (n_threads > 1) {
    if (in.is_bgzf()) tp.set_io(in);
    tp.set_io(out);
  }

  static constexpr uint32_t output_buffer_size = 1024;
  vector<char> buf(output_buffer_size, '\0');

  kstring_t line{0, 0, nullptr};
  const int ret = ks_resize(&line, output_buffer_size);
  if (ret) throw runtime_error("failed to acquire buffer");

  string chrom_name;
  uint32_t nm_sz{};
  int32_t prev_chr_id = -1;
  uint64_t pos = num_lim<uint64_t>::max();

  // ADS: this is probably a poor strategy since we already would know
  // the index of the chrom sequence in the vector.
  chrom_itr_t ch_itr;

  while (getline(in, line)) {
    if (is_counts_header_line(line.s)) {
      process_header_line(name_to_id, chrom_sizes, line, out);
      continue;  // ADS: early loop exit
    }

    if (!std::isdigit(line.s[0])) {  // check if we have a chrom line

      if (!require_covered && pos != num_lim<uint64_t>::max())
        write_missing_cpg(nm_sz, *ch_itr, pos + 1, size(*ch_itr), buf, out);

      chrom_name = string{line.s};
      nm_sz = size(chrom_name);
      const int32_t chr_id = get_chrom_id(name_to_id, chrom_name);

      if (add_missing_chroms)
        write_all_cpgs(verbose, prev_chr_id, chr_id, names, chroms, buf, out);

      ch_itr = get_chrom(chrom_lookup, chrom_name);
      pos = 0;
      prev_chr_id = chr_id;
      if (verbose) cerr << "processing: " << chrom_name << endl;

      auto res = copy(cbegin(chrom_name), cend(chrom_name), buf.data());
      *res = '\t';
    }
    else {
      uint32_t pos_step = 0, n_meth = 0, n_unmeth = 0;
      const auto end_line = line.s + line.l;
      auto res = from_chars(line.s, end_line, pos_step);
      res = from_chars(res.ptr + 1, end_line, n_meth);
      res = from_chars(res.ptr + 1, end_line, n_unmeth);

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

    string outfile;
    string chroms_file;
    const string description =
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
    std::vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl
           << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    if (leftover_args.size() != 1) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    if (require_covered && add_missing_chroms) {
      cerr << "options mutually exclusive: reads and missing" << endl;
      return EXIT_FAILURE;
    }
    const string filename(leftover_args.front());
    /****************** END COMMAND LINE OPTIONS *****************/

    if (require_covered && add_missing_chroms) {
      cerr << "options mutually exclusive: reads and missing" << endl;
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
  catch (const std::runtime_error &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
