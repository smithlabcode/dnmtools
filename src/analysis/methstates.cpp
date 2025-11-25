/* methstates: a program for converting read sequences in SAM format
 * files into methylation states at CpGs covered by those reads
 *
 * Copyright (C) 2011-2022 University of Southern California and
 *                         Andrew D. Smith
 *
 * Authors: Andrew D. Smith and Masaru Nakajima
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

#include "OptionParser.hpp"
#include "bam_record_utils.hpp"
#include "dnmt_error.hpp"
#include "smithlab_os.hpp"

#include <bamxx.hpp>

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <iterator>
#include <limits>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

using std::cerr;
using std::cout;
using std::lower_bound;
using std::string;
using std::unordered_map;
using std::unordered_set;
using std::vector;

using bamxx::bam_rec;

static const char b2c[] = "TNGNNNCNNNNNNNNNNNNA";  // NOLINT(*-avoid-c-arrays)

struct quick_buf : public std::ostringstream,
                   public std::basic_stringbuf<char> {
  // ADS: consider putting this class in a header somewhere
  quick_buf() { static_cast<std::basic_ios<char> &>(*this).rdbuf(this); }

  void
  clear() {
    setp(pbase(), pbase());
  }

  char const *
  c_str() {
    *pptr() = '\0';
    return pbase();
  }
};

template <class BidirIt, class OutputIt>
// constexpr // since C++20
OutputIt
revcomp_copy(BidirIt first, BidirIt last, OutputIt d_first) {
  for (; first != last; ++d_first)
    *d_first = b2c[*(--last) - 'A'];  // NOLINT(*-constant-array-index)
  return d_first;
}

inline static bool
is_cpg(const string &s, const uint64_t idx) {
  return s[idx] == 'C' && s[idx + 1] == 'G';
}

static void
collect_cpgs(const string &s, vector<uint64_t> &cpgs) {
  cpgs.clear();
  const uint64_t lim = std::size(s) - 1;
  for (auto i = 0u; i < lim; ++i)
    if (is_cpg(s, i))
      cpgs.push_back(i);
}

static bool
convert_meth_states_pos(const vector<uint64_t> &cpgs,
                        const bamxx::bam_header &hdr, const bam_rec &aln,
                        uint64_t &first_cpg_index, string &states) {
  states.clear();

  const uint64_t seq_start = get_pos(aln);
  const uint64_t width = rlen_from_cigar(aln);
  const uint64_t seq_end = seq_start + width;

  string seq_str;
  get_seq_str(aln, seq_str);
  apply_cigar(aln, seq_str, 'N');

  if (std::size(seq_str) != width)
    throw dnmt_error("bad sam record format: " + to_string(hdr, aln));

  // get the first cpg site equal to or large than seq_start
  auto cpg_itr = lower_bound(begin(cpgs), end(cpgs), seq_start);
  auto first_cpg_itr = end(cpgs);

  if (cpg_itr == end(cpgs))
    return false;

  for (; cpg_itr != end(cpgs) && *cpg_itr < seq_end; cpg_itr++) {
    const char x = seq_str[*cpg_itr - seq_start];
    states += (x == 'T') ? 'T' : ((x == 'C') ? 'C' : 'N');
    if (first_cpg_itr == end(cpgs))
      first_cpg_itr = cpg_itr;
  }

  if (first_cpg_itr != end(cpgs))
    first_cpg_index = distance(begin(cpgs), first_cpg_itr);

  return states.find_first_of("CT") != string::npos;
}

static bool
convert_meth_states_neg(const vector<uint64_t> &cpgs,
                        const bamxx::bam_header &hdr, const bam_rec &aln,
                        uint64_t &first_cpg_index, string &states) {
  /* ADS: the "revcomp" on the read sequence is needed for the cigar
     to be applied, since the cigar is relative to the genome
     coordinates and not the read's sequence. But the read sequence
     may is assumed to have been T-rich to begin with, so it becomes
     A-rich. And the position of the C in the CpG becomes the G
     position.
   */

  states.clear();

  const uint64_t seq_start = get_pos(aln);
  const uint64_t width = rlen_from_cigar(aln);
  const uint64_t seq_end = seq_start + width;

  string orig_seq;
  get_seq_str(aln, orig_seq);

  string seq_str;
  seq_str.resize(orig_seq.size());
  revcomp_copy(begin(orig_seq), end(orig_seq), begin(seq_str));
  apply_cigar(aln, seq_str, 'N');

  if (seq_str.size() != width)
    throw dnmt_error("bad sam record format: " + to_string(hdr, aln));

  // get the first cpg site equal to or large than seq_start - 1
  // the -1 is because we look for G in the read corresponding to a
  // CpG in chromosome, which are indexed in cpgs based on the position of C
  auto cpg_itr =
    lower_bound(begin(cpgs), end(cpgs), seq_start > 0 ? seq_start - 1 : 0);
  auto first_cpg_itr = end(cpgs);

  if (cpg_itr == end(cpgs)) {
    return false;
  }
  else {
    for (; cpg_itr != end(cpgs) && *cpg_itr < seq_end - 1; cpg_itr++) {
      const char x = seq_str[*cpg_itr - seq_start + 1];
      states += (x == 'G') ? 'C' : ((x == 'A') ? 'T' : 'N');
      if (first_cpg_itr == end(cpgs))
        first_cpg_itr = cpg_itr;
    }
  }

  if (first_cpg_itr != end(cpgs)) {
    first_cpg_index = distance(begin(cpgs), first_cpg_itr);
  }

  return states.find_first_of("CT") != string::npos;
}

static void
get_chrom(const string &chrom_name, const vector<string> &all_chroms,
          const unordered_map<string, uint64_t> &chrom_lookup, string &chrom) {
  auto the_chrom = chrom_lookup.find(chrom_name);
  if (the_chrom == end(chrom_lookup))
    throw dnmt_error("could not find chrom: " + chrom_name);

  chrom = all_chroms[the_chrom->second];
  if (chrom.empty())
    throw dnmt_error("problem with chrom: " + chrom_name);
}

int
main_methstates(int argc, char *argv[]) {  // NOLINT(*-avoid-c-arrays)
  try {
    const auto description =
      R"(Convert mapped reads in SAM format into a format that indicates binary
      sequences of methylation states in each read, indexed by the identity
      of the CpG they cover, along with the chromosome. Only reads that
      cover a CpG site are included in the output. All output is relative to
      the positive reference strand. This format is used as input to other
      tools, and is not intended to be human-interpretable. All chromosome
      sequences are loaded at once.)";

    bool VERBOSE = false;
    bool compress_output = false;

    string chrom_file;
    string outfile("-");
    int n_threads = 1;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(argv[0],  // NOLINT(*-pointer-arithmetic)
                           description, "<sam-file>");
    opt_parse.add_opt("output", 'o', "output file name", false, outfile);
    opt_parse.add_opt("chrom", 'c', "fasta format reference genome file", true,
                      chrom_file);
    opt_parse.add_opt("threads", 't', "threads to use for reading input", false,
                      n_threads);
    opt_parse.add_opt("zip", 'z', "output gzip format", false, compress_output);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);

    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << '\n'
           << opt_parse.about_message() << '\n';
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      cerr << opt_parse.about_message() << '\n';
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << '\n';
      return EXIT_SUCCESS;
    }
    if (leftover_args.size() != 1) {
      cerr << opt_parse.help_message() << '\n';
      return EXIT_SUCCESS;
    }
    const string mapped_reads_file = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/

    if (n_threads < 0)
      throw dnmt_error("thread count cannot be negative");

    /* first load in all the chromosome sequences and names, and make
       a map from chromosome name to the location of the chromosome
       itself */
    vector<string> all_chroms, chrom_names;
    read_fasta_file_short_names(chrom_file, chrom_names, all_chroms);
    for (auto &&i : all_chroms)
      transform(begin(i), end(i), begin(i),
                [](const char c) { return std::toupper(c); });

    unordered_map<string, uint64_t> chrom_lookup;
    for (uint64_t i = 0; i < chrom_names.size(); ++i)
      chrom_lookup[chrom_names[i]] = i;

    if (VERBOSE)
      cerr << "n_chroms: " << all_chroms.size() << '\n';

    bamxx::bam_tpool tp(n_threads);  // declared first; destroyed last

    bamxx::bam_in in(mapped_reads_file);
    if (!in)
      throw dnmt_error("cannot open input file " + mapped_reads_file);
    bamxx::bam_header hdr(in);
    if (!hdr)
      throw dnmt_error("cannot read heade" + mapped_reads_file);

    // open the output file
    const string output_mode = compress_output ? "w" : "wu";
    bamxx::bgzf_file out(outfile, output_mode);
    if (!out)
      throw dnmt_error("error opening output file: " + outfile);

    /* set the threads for the input file decompression */
    if (n_threads > 1) {
      tp.set_io(in);
      tp.set_io(out);
    }

    vector<uint64_t> cpgs;

    unordered_set<string> chroms_seen;
    string chrom_name;
    string chrom;

    quick_buf buf;

    // iterate over records/reads in the SAM file, sequentially
    // processing each before considering the next
    bam_rec aln;
    while (in.read(hdr, aln)) {
      // get the correct chrom if it has changed
      if (string(sam_hdr_tid2name(hdr, aln)) != chrom_name) {
        chrom_name = sam_hdr_tid2name(hdr, aln);

        // make sure all reads from same chrom are contiguous in the file
        if (chroms_seen.find(chrom_name) != end(chroms_seen))
          throw dnmt_error("chroms out of order (check SAM file sorted)");

        if (VERBOSE)
          cerr << "processing " << chrom_name << '\n';

        get_chrom(chrom_name, all_chroms, chrom_lookup, chrom);
        collect_cpgs(chrom, cpgs);
      }

      uint64_t first_cpg_index = std::numeric_limits<uint64_t>::max();
      string seq;

      const bool has_cpgs =
        bam_is_rev(aln)
          ? convert_meth_states_neg(cpgs, hdr, aln, first_cpg_index, seq)
          : convert_meth_states_pos(cpgs, hdr, aln, first_cpg_index, seq);

      if (has_cpgs) {
        buf.clear();
        buf << sam_hdr_tid2name(hdr, aln) << '\t' << first_cpg_index << '\t'
            << seq << '\n';
        if (!out.write(buf.c_str(), buf.tellp())) {
          cerr << "failure writing output" << '\n';
          return EXIT_FAILURE;
        }
      }
    }
  }
  catch (const std::exception &e) {
    cerr << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
