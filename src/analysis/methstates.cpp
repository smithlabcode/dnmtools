/* methstates: a program for converting read sequences in SAM format
 * files into methylation states at CpGs covered by those reads
 *
 * Copyright (C) 2011-2022 University of Southern California and
 *                         Andrew D. Smith
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

#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <unordered_set>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "htslib_wrapper.hpp"
#include "sam_record.hpp"
#include "cigar_utils.hpp"

#include "bam_record_utils.hpp"

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::unordered_map;
using std::unordered_set;
using std::runtime_error;
using std::lower_bound;


static const char b2c[] = "TNGNNNCNNNNNNNNNNNNA";

template<class BidirIt, class OutputIt>
// constexpr // since C++20
OutputIt revcomp_copy(BidirIt first, BidirIt last, OutputIt d_first) {
  for (; first != last; ++d_first)
    *d_first = b2c[*(--last) - 'A'];
  return d_first;
}


// template<class BidirIt>
// constexpr // since C++20
// void revcomp_into(BidirIt first, BidirIt last)
// {
//     using iter_cat = typename std::iterator_traits<BidirIt>::iterator_category;

//     // Tag dispatch, e.g. calling reverse_impl(first, last, iter_cat()),
//     // can be used in C++14 and earlier modes.
//     if constexpr (std::is_base_of_v<std::random_access_iterator_tag, iter_cat>)
//     {
//         if (first == last)
//             return;

//         for (--last; first < last; (void)++first, --last)
//             std::iter_swap(first, last);
//     }
//     else
//         while (first != last && first != --last)
//             std::iter_swap(first++, last);
// }


inline static bool
is_cpg(const string &s, const size_t idx) {
  return s[idx] == 'C' && s[idx + 1] == 'G';
}

static void
collect_cpgs(const string &s, unordered_map<size_t, size_t> &cpgs) {
  cpgs.clear();
  const size_t lim = s.length() - 1;
  size_t cpg_count = 0;
  for (size_t i = 0; i < lim; ++i)
    if (is_cpg(s, i))
      cpgs[i] = cpg_count++;
}

static void
collect_cpgs(const string &s, 
             vector<size_t> &cpgs) {
  cpgs.clear();
  const size_t lim = s.length() - 1;
  for (size_t i = 0; i < lim; ++i) {
    if (is_cpg(s, i)) 
      cpgs.push_back(i);
  }
}

// void
// local_apply_cigar(const string &cigar, string &to_inflate,
//             const char inflation_symbol) {
//   std::istringstream iss(cigar);

//   string inflated_seq;
//   size_t n;
//   char op;
//   size_t i = 0;
//   auto to_inflate_beg = std::begin(to_inflate);
//   while (iss >> n >> op) {
//     if (consumes_reference(op) && consumes_query(op)) {
//       inflated_seq.append(to_inflate_beg + i, to_inflate_beg + i + n);
//       i += n;
//     }
//     else if (consumes_query(op)) {
//       // no addition of symbols to query
//       i += n;
//     }
//     else if (consumes_reference(op)) {
//       inflated_seq.append(n, inflation_symbol);
//       // no increment of index within query
//     }
//   }

//   // sum of total M/I/S/=/X/N operations must equal length of seq
//   const size_t orig_len = to_inflate.length();
//   if (i != orig_len)
//     throw runtime_error("inconsistent number of qseq ops in cigar: " +
//                         to_inflate + " "  + cigar + " " +
//                         to_string(i) + " " +
//                         to_string(orig_len));
//   to_inflate.swap(inflated_seq);
// }

static bool
convert_meth_states_pos(const string &chrom,
                        const vector<size_t> &cpgs,
                        const bam_header &hdr, 
                        const bam_rec &aln,
                        size_t &first_cpg_index, string &states) {

  states.clear();

  const size_t seq_start = get_pos(aln);  
  const size_t width = rlen_from_cigar(aln);
  const size_t seq_end = seq_start + width;   

  string seq_str;
  get_seq_str(aln, seq_str);
  apply_cigar(aln, seq_str, 'N');
  if (seq_str.size() != width) {
    throw runtime_error("bad sam record format: " + to_string(hdr, aln));
  }

  // get the first cpg site equal to or large than seq_start
  auto cpg_itr = lower_bound(begin(cpgs), end(cpgs), seq_start);
  auto first_cpg_itr = end(cpgs);

  if (cpg_itr == end(cpgs)) {
    return false;
  } else {
    for (; cpg_itr != end(cpgs) && *cpg_itr < seq_end; cpg_itr++) {
      const char x = seq_str[*cpg_itr - seq_start];
      states += (x == 'C') ? 'C' : ((x == 'T') ? 'T' : 'N');
      if (first_cpg_itr == end(cpgs))
        first_cpg_itr = cpg_itr;
    }
  }

  if (first_cpg_itr != end(cpgs)) {
    first_cpg_index = distance(begin(cpgs), first_cpg_itr);
  }

  return states.find_first_of("CT") != string::npos;
}


static bool
convert_meth_states_pos(const string &chrom,
                        const unordered_map<size_t, size_t> &cpgs,
                        const bam_header &hdr, 
                        const bam_rec &aln,
                        size_t &start_pos, string &states) {

  states.clear();

  const size_t width = rlen_from_cigar(aln);
  const size_t offset = get_pos(aln);  

  string seq_str;
  get_seq_str(aln, seq_str);
  apply_cigar(aln, seq_str, 'N');
  if (seq_str.size() != width) {
    throw runtime_error("bad sam record format: " + to_string(hdr, aln));
  }

  const auto beg_chrom = begin(chrom);
  const auto end_chrom = end(chrom);
  auto chrom_itr = beg_chrom + offset;
  auto first_cpg = end_chrom;

  // ADS: below the "-1" in the std::min is to ensure dinuc exists
  const auto beg_seq = begin(seq_str);
  const auto seq_lim = beg_seq + std::min(width, chrom.length() - 1 - offset);

  for (auto seq_itr = beg_seq; seq_itr != seq_lim; ++chrom_itr, ++seq_itr) {
    if (*chrom_itr == 'C' && *(chrom_itr + 1) == 'G') {
      const char x = *seq_itr;
      states += (x == 'C') ? 'C' : ((x == 'T') ? 'T' : 'N');
      if (first_cpg == end_chrom)
        first_cpg = chrom_itr;
    }
  }

  if (first_cpg != end_chrom) {
    const auto the_cpg = cpgs.find(distance(beg_chrom, first_cpg));
    if (the_cpg == end(cpgs))
      throw runtime_error("cannot locate site on positive strand: " +
                          to_string(hdr, aln));
    start_pos = the_cpg->second;
  }
  return states.find_first_of("CT") != string::npos;
}


static bool
convert_meth_states_neg(const string &chrom,
                        const vector<size_t> &cpgs,
                        const bam_header &hdr,
                        const bam_rec &aln,
                        size_t &first_cpg_index, string &states) {
  /* ADS: the "revcomp" on the read sequence is needed for the cigar
     to be applied, since the cigar is relative to the genome
     coordinates and not the read's sequence. But the read sequence
     may is assumed to have been T-rich to begin with, so it becomes
     A-rich. And the position of the C in the CpG becomes the G
     position.
   */

  states.clear();

  const size_t seq_start = get_pos(aln);  
  const size_t width = rlen_from_cigar(aln);
  const size_t seq_end = seq_start + width;   

  string orig_seq;
  get_seq_str(aln, orig_seq);
  string seq_str;
  seq_str.resize(orig_seq.size());
  revcomp_copy(begin(orig_seq), end(orig_seq), begin(seq_str));
  apply_cigar(aln, seq_str, 'N');
  if (seq_str.size() != width) {
    throw runtime_error("bad sam record format: " + to_string(hdr, aln));
  }

  // get the first cpg site equal to or large than seq_start - 1
  // the -1 is because we look for G in the read corresponding to a 
  // CpG in chromosome, which are indexed in cpgs based on the position of C 
  auto cpg_itr = lower_bound(begin(cpgs), end(cpgs), 
                             seq_start > 0 ? seq_start - 1 : 0);
  auto first_cpg_itr = end(cpgs);

  if (cpg_itr == end(cpgs)) {
    return false;
  } else {
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

static bool
convert_meth_states_neg(const string &chrom,
                        const unordered_map<size_t, size_t> &cpgs,
                        const bam_header &hdr,
                        const bam_rec &aln,
                        size_t &start_pos, string &states) {
  /* ADS: the "revcomp" on the read sequence is needed for the cigar
     to be applied, since the cigar is relative to the genome
     coordinates and not the read's sequence. But the read sequence
     may is assumed to have been T-rich to begin with, so it becomes
     A-rich. And the position of the C in the CpG becomes the G
     position.
   */

  states.clear();

  const size_t width = rlen_from_cigar(aln);
  const size_t offset = get_pos(aln);  

  string orig_seq;
  get_seq_str(aln, orig_seq);
  string the_seq;
  the_seq.resize(orig_seq.size());
  revcomp_copy(begin(orig_seq), end(orig_seq), begin(the_seq));
  apply_cigar(aln, the_seq, 'N');
  if (the_seq.size() != width) {
    throw runtime_error("bad sam record format: " + to_string(hdr, aln));
  }

  const auto beg_chrom = begin(chrom);
  const auto end_chrom = end(chrom);
  auto chrom_itr = beg_chrom + ((offset > 0) ? offset - 1 : 0);
  auto first_cpg = end_chrom;

  // ADS: should apply_cigar make the std::min below irrelevant?
  const auto beg_seq = begin(the_seq);
  const auto seq_lim = beg_seq + std::min(width, chrom.length() - offset);
  auto seq_itr = beg_seq + ((offset > 0) ? 0 : 1);

  for (; seq_itr != seq_lim; ++chrom_itr, ++seq_itr) {
    if (*chrom_itr == 'C' && *(chrom_itr + 1) == 'G') {
      const char x = *seq_itr;
      states += (x == 'G') ? 'C' : ((x == 'A') ? 'T' : 'N');
      if (first_cpg == end_chrom)
        first_cpg = chrom_itr;
    }
  }

  if (first_cpg != end_chrom) {
    const auto the_cpg = cpgs.find(distance(beg_chrom, first_cpg));
    if (the_cpg == end(cpgs))
      throw runtime_error("cannot locate site on negative strand: " +
                          to_string(hdr, aln));
    start_pos = the_cpg->second;
  }
  return states.find_first_of("CT") != string::npos;
}


static void
get_chrom(const string &chrom_name,
          const vector<string> &all_chroms,
          const unordered_map<string, size_t> &chrom_lookup,
          string &chrom) {

  auto the_chrom = chrom_lookup.find(chrom_name);
  if (the_chrom == end(chrom_lookup))
    throw runtime_error("could not find chrom: " + chrom_name);

  chrom = all_chroms[the_chrom->second];
  if (chrom.empty())
    throw runtime_error("problem with chrom: " + chrom_name);
}


int
main_methstates(int argc, const char **argv) {

  try {

    const string description =
      "Convert mapped reads in SAM format into a format that indicates binary \
      sequences of methylation states in each read, indexed by the identity   \
      of the CpG they cover, along with the chromosome. Only reads that       \
      cover a CpG site are included in the output. All output is relative to  \
      the positive reference strand. This format is used as input to other    \
      tools, and is not intended to be human-interpretable. All chromosome    \
      sequences are loaded at once.";

    bool VERBOSE = false;

    string chrom_file;
    string outfile;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(argv[0], description, "<sam-file>");
    opt_parse.add_opt("output", 'o', "output file name", false, outfile);
    opt_parse.add_opt("chrom", 'c', "fasta format reference genome file",
                      true , chrom_file);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);

    vector<string> leftover_args;
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
    const string mapped_reads_file = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/

    /* first load in all the chromosome sequences and names, and make
       a map from chromosome name to the location of the chromosome
       itself */
    vector<string> all_chroms, chrom_names;
    read_fasta_file_short_names(chrom_file, chrom_names, all_chroms);
    for (auto &&i: all_chroms)
      transform(begin(i), end(i), begin(i),
                [](const char c) { return std::toupper(c); });

    unordered_map<string, size_t> chrom_lookup;
    for (size_t i = 0; i < chrom_names.size(); ++i)
      chrom_lookup[chrom_names[i]] = i;

    if (VERBOSE)
      cerr << "n_chroms: " << all_chroms.size() << endl;

    bam_infile in(mapped_reads_file);
    if (!in)
      throw runtime_error("cannot open input file " + mapped_reads_file);
    bam_header hdr(in);
    if (!hdr)
      throw runtime_error("cannot read heade" + mapped_reads_file);

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());

    // for the current chrom, this maps cpg index to cpg position in
    // the chrom
    //unordered_map<size_t, size_t> cpgs;
    vector<size_t> cpgs;

    unordered_set<string> chroms_seen;
    string chrom_name;
    string chrom;

    // iterate over records/reads in the SAM file, sequentially
    // processing each before considering the next
    bam_rec aln;
    while (in.read(hdr, aln)) {

      // get the correct chrom if it has changed
      if (string(sam_hdr_tid2name(hdr, aln)) != chrom_name) {
        chrom_name = sam_hdr_tid2name(hdr, aln);

        // make sure all reads from same chrom are contiguous in the file
        if (chroms_seen.find(chrom_name) != end(chroms_seen))
          throw runtime_error("chroms out of order (check SAM file sorted)");

        if (VERBOSE)
          cerr << "processing " << chrom_name << endl;

        get_chrom(chrom_name, all_chroms, chrom_lookup, chrom);
        collect_cpgs(chrom, cpgs);
      }

      size_t first_cpg_index = std::numeric_limits<size_t>::max();
      string seq;

      const bool has_cpgs = bam_is_rev(aln) ?
        convert_meth_states_neg(chrom, cpgs, hdr, aln, first_cpg_index, seq) :
        convert_meth_states_pos(chrom, cpgs, hdr, aln, first_cpg_index, seq);

      if (has_cpgs)
        out << sam_hdr_tid2name(hdr, aln) << '\t'
            << first_cpg_index << '\t'
            << seq << '\n';
    }
  }
  catch (const std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
