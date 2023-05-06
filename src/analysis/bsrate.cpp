/* bsrate: a program for determining the rate of bisulfite conversion
 * in a bisulfite sequencing experiment
 *
 * Copyright (C) 2014-2023 University of Southern California and
 *                         Andrew D. Smith
 *
 * Authors: Andrew D. Smith and Guilherme Sena
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
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"
#include "zlib_wrapper.hpp"
#include "bsutils.hpp"
#include "cigar_utils.hpp"
#include "htslib_wrapper.hpp"
#include "sam_record.hpp"

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <unordered_set>

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::max;
using std::accumulate;
using std::to_string;
using std::unordered_map;
using std::unordered_set;
using std::runtime_error;

inline bool
is_rc(const sam_rec &aln) {
  return check_flag(aln, samflags::read_rc);
}

static void
flip_strand(sam_rec &aln) {
  if (check_flag(aln, samflags::read_rc))
    unset_flag(aln, samflags::read_rc);
  else
    set_flag(aln, samflags::read_rc);
}

static void
revcomp(sam_rec &aln) {
  flip_strand(aln); // set strand to opposite of current value
  revcomp_inplace(aln.seq); // reverse complement sequence
  std::reverse(begin(aln.qual), end(aln.qual)); // and quality scores
}


static void
count_states_pos(const bool INCLUDE_CPGS, const string &chrom,
                 const sam_rec &aln,
                 vector<size_t> &unconv, vector<size_t> &conv,
                 vector<size_t> &err, size_t &hanging) {

  const size_t width = cigar_rseq_ops(aln.cigar);
  const size_t offset = aln.pos - 1;

  string seq(aln.seq);
  apply_cigar(aln.cigar, seq);

  assert(seq.size() == width);

  size_t position = offset;
  if (chrom.length() < offset) // at least one bp of read on chr
    throw runtime_error("read mapped off chrom:\n" +
                        to_string(chrom.length()));
  for (size_t i = 0; i < width; ++i, ++position) {
    if (position >= chrom.length()) // some overhang
      ++hanging;

    if (is_cytosine(chrom[position]) &&
        (!is_guanine(chrom[position+1]) ||
         position == chrom.length() ||
         INCLUDE_CPGS)) {

      if (is_cytosine(seq[i])) ++unconv[i];
      else if (is_thymine(seq[i])) ++conv[i];
      else if (toupper(seq[i]) != 'N')
        ++err[i];
    }
  }
}


static void
count_states_neg(const bool INCLUDE_CPGS, const string &chrom,
                 const sam_rec &aln,
                 vector<size_t> &unconv, vector<size_t> &conv,
                 vector<size_t> &err, size_t &hanging) {

  const size_t width = cigar_rseq_ops(aln.cigar);
  const size_t offset = aln.pos - 1;

  // ADS: reverse complement twice because the cigar is applied
  // starting relative to the reference.
  string seq(aln.seq);
  revcomp_inplace(seq);
  apply_cigar(aln.cigar, seq);
  revcomp_inplace(seq);

  assert(seq.size() == width);

  size_t position = offset + width - 1;
  assert(offset < chrom.length()); // at least one bp of read on chr
  for (size_t i = 0; i < width; ++i, --position) {
    if (position >= chrom.length())
      ++hanging;
    if (is_guanine(chrom[position]) &&
        (position == 0 ||
         !is_cytosine(chrom[position-1]) ||
         INCLUDE_CPGS)) {

      if (is_cytosine(seq[i])) ++unconv[i];
      else if (is_thymine(seq[i])) ++conv[i];
      else if (toupper(seq[i]) != 'N')
        ++err[i];
    }
  }
}


static void
write_output(const string &outfile,
             const vector<size_t> &ucvt_count_p,
             const vector<size_t> &cvt_count_p,
             const vector<size_t> &ucvt_count_n,
             const vector<size_t> &cvt_count_n,
             const vector<size_t> &err_p, const vector<size_t> &err_n) {

  // Get some totals first
  const size_t pos_cvt = accumulate(begin(cvt_count_p), end(cvt_count_p), 0ul);
  const size_t neg_cvt = accumulate(begin(cvt_count_n), end(cvt_count_n), 0ul);
  const size_t total_cvt = pos_cvt + neg_cvt;

  const size_t pos_ucvt =
    accumulate(begin(ucvt_count_p), end(ucvt_count_p), 0ul);
  const size_t neg_ucvt =
    accumulate(begin(ucvt_count_n), end(ucvt_count_n), 0ul);
  const size_t total_ucvt = pos_ucvt + neg_ucvt;

  std::ofstream of;
  if (!outfile.empty()) of.open(outfile.c_str());
  std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

  out << "OVERALL CONVERSION RATE = "
      << static_cast<double>(total_cvt)/(total_cvt + total_ucvt) << endl
      << "POS CONVERSION RATE = "
      << static_cast<double>(pos_cvt)/(pos_cvt + pos_ucvt) << '\t'
      << std::fixed << static_cast<size_t>(pos_cvt + pos_ucvt) << endl
      << "NEG CONVERSION RATE = "
      << static_cast<double>(neg_cvt)/(neg_cvt + neg_ucvt) << '\t'
      << std::fixed << static_cast<size_t>(neg_cvt + neg_ucvt) << endl;

  out << "BASE" << '\t'
      << "PTOT" << '\t'
      << "PCONV" << '\t'
      << "PRATE" << '\t'
      << "NTOT" << '\t'
      << "NCONV" << '\t'
      << "NRATE" << '\t'
      << "BTHTOT" << '\t'
      << "BTHCONV" << '\t'
      << "BTHRATE" << '\t'
      << "ERR" << '\t'
      << "ALL" << '\t'
      << "ERRRATE"  << endl;

  // Figure out how many positions to print in the output, capped at 1000
  size_t output_len =
    (ucvt_count_p.size() > 1000) ? 1000 : ucvt_count_p.size();

  while (output_len > 0 &&
         (ucvt_count_p[output_len-1] + cvt_count_p[output_len-1] +
          ucvt_count_n[output_len-1] + cvt_count_n[output_len-1] == 0))
    --output_len;

  // Now actually output the results
  static const size_t precision_val = 5;
  for (size_t i = 0; i < output_len; ++i) {

    const size_t total_p = ucvt_count_p[i] + cvt_count_p[i];
    const size_t total_n = ucvt_count_n[i] + cvt_count_n[i];
    const size_t total_valid = total_p + total_n;
    out << (i + 1) << "\t";

    out.precision(precision_val);
    out << total_p << '\t' << cvt_count_p[i] << '\t'
        << static_cast<double>(cvt_count_p[i])/max(size_t(1ul), total_p)
        << '\t';

    out.precision(precision_val);
    out << total_n << '\t' << cvt_count_n[i] << '\t'
        << static_cast<double>(cvt_count_n[i])/max(size_t(1ul), total_n)
        << '\t';

    const double total_cvt = cvt_count_p[i] + cvt_count_n[i];
    out.precision(precision_val);
    out << static_cast<size_t>(total_valid)
        << '\t' << cvt_count_p[i] + cvt_count_n[i] << '\t'
        << total_cvt/max(1ul, total_valid) << '\t';

    const double total_err = err_p[i] + err_n[i];
    out.precision(precision_val);
    const size_t total = total_valid + err_p[i] + err_n[i];
    out << err_p[i] + err_n[i] << '\t' << static_cast<size_t>(total) << '\t'
        << total_err/max(1ul, total) << endl;
  }
}


int
main_bsrate(int argc, const char **argv) {

  try {

    // ASSUMED MAXIMUM LENGTH OF A FRAGMENT
    static const size_t output_size = 10000;

    bool VERBOSE = false;
    bool INCLUDE_CPGS = false;
    bool reads_are_a_rich = false;

    string chroms_file;
    string outfile;
    string seq_to_use; // use only this chrom/sequence in the analysis


    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "Program to compute the "
                           "BS conversion rate from BS-seq "
                           "reads mapped to a genome",
                           "-c <chroms> <mapped-reads>");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)",
                      false, outfile);
    opt_parse.add_opt("chrom", 'c', "File of chromosome sequences (FASTA)",
                      true, chroms_file);
    opt_parse.add_opt("all", 'N', "count all Cs (including CpGs)",
                      false , INCLUDE_CPGS);
    opt_parse.add_opt("seq", '\0', "use only this sequence (e.g. chrM)",
                      false , seq_to_use);
    opt_parse.add_opt("a-rich", 'A', "reads are A-rich",
                      false, reads_are_a_rich);
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

    vector<string> chroms;
    vector<string> names;
    read_fasta_file_short_names(chroms_file, names, chroms);
    unordered_map<string, size_t> chrom_lookup;
    for (size_t i = 0; i < chroms.size(); ++i)
      chrom_lookup[names[i]] = i;

    if (VERBOSE)
      cerr << "[n chroms in reference: " << chroms.size() << "]" << endl;

    igzfstream in(mapped_reads_file);
    if (!in)
      throw runtime_error("cannot open file: " + mapped_reads_file);

    vector<size_t> unconv_count_pos(output_size, 0ul);
    vector<size_t> conv_count_pos(output_size, 0ul);
    vector<size_t> unconv_count_neg(output_size, 0ul);
    vector<size_t> conv_count_neg(output_size, 0ul);
    vector<size_t> err_pos(output_size, 0ul);
    vector<size_t> err_neg(output_size, 0ul);

    size_t chrom_idx = 0;
    string chrom_name;
    size_t hanging = 0;

    bool use_this_chrom = seq_to_use.empty();

    SAMReader sam_reader(mapped_reads_file);
    sam_rec aln;

    unordered_set<string> chroms_seen;
    while (sam_reader >> aln) {

      if (reads_are_a_rich)
        revcomp(aln);

      // get the correct chrom if it has changed
      if (aln.rname != chrom_name) {
        // make sure all reads from same chrom are contiguous in the file
        if (chroms_seen.find(aln.rname) != end(chroms_seen))
          throw runtime_error("chroms out of order in mapped reads file");

        if (VERBOSE)
          cerr << "processing chromosome " << aln.rname << "\n";

        chroms_seen.insert(aln.rname);

        auto chrom_itr = chrom_lookup.find(aln.rname);
        if (chrom_itr == end(chrom_lookup))
          throw runtime_error("could not find chrom: " + aln.rname);

        chrom_idx = chrom_itr->second;
        chrom_name = aln.rname;
        use_this_chrom = seq_to_use.empty() || chrom_name == seq_to_use;
      }

      if (use_this_chrom) {
        // do the work for this mapped read
        if (is_rc(aln))
          count_states_neg(INCLUDE_CPGS, chroms[chrom_idx], aln,
                           unconv_count_neg, conv_count_neg, err_neg, hanging);
        else
          count_states_pos(INCLUDE_CPGS, chroms[chrom_idx], aln,
                           unconv_count_pos, conv_count_pos, err_pos, hanging);
      }
    }
    write_output(outfile, unconv_count_pos,
                 conv_count_pos, unconv_count_neg,
                 conv_count_neg, err_pos, err_neg);

    if (hanging > 0) // some overhanging reads
      cerr << "Warning: hanging reads detected at chrom ends "
           << "(N=" << hanging<< ")" << endl
           << "High numbers of hanging reads suggest mismatch "
           << "between assembly provided here and that used for mapping"
           << endl;
  }
  catch (const runtime_error &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
