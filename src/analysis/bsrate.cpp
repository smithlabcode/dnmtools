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

#include <algorithm>
#include <fstream>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

#include "OptionParser.hpp"
#include "bam_record_utils.hpp"
#include "bsutils.hpp"
#include "dnmt_error.hpp"
#include "smithlab_utils.hpp"

#include <bamxx.hpp>

using std::accumulate;
using std::cerr;
using std::cout;
using std::endl;
using std::max;
using std::numeric_limits;
using std::runtime_error;
using std::string;
using std::unordered_map;
using std::unordered_set;
using std::vector;

using bamxx::bam_rec;

static void
count_states_pos(const bool INCLUDE_CPGS, const string &chrom,
                 const bam_rec &aln, vector<size_t> &unconv,
                 vector<size_t> &conv, vector<size_t> &err, size_t &hanging) {
  /* iterate through reference, query/read and fragment */
  const auto seq = bam_get_seq(aln);
  const auto beg_cig = bam_get_cigar(aln);
  const auto end_cig = beg_cig + get_n_cigar(aln);
  auto rpos = get_pos(aln);
  auto qpos = 0;
  auto fpos = 0;

  const decltype(rpos) chrom_lim = chrom.size() - 1;

  for (auto c_itr = beg_cig; c_itr != end_cig; ++c_itr) {
    const auto op = bam_cigar_op(*c_itr);
    const auto n = bam_cigar_oplen(*c_itr);
    if (cigar_eats_ref(op) && cigar_eats_query(op)) {
      const decltype(qpos) end_qpos = qpos + n;
      for (; qpos < end_qpos; ++qpos, ++rpos, ++fpos) {
        if (rpos > chrom_lim) ++hanging;
        if (is_cytosine(chrom[rpos]) &&
            (rpos >= chrom_lim || !is_guanine(chrom[rpos + 1]) ||
             INCLUDE_CPGS)) {
          const auto qc = seq_nt16_str[bam_seqi(seq, qpos)];
          if (qc == 'C')
            ++unconv[fpos];
          else if (qc == 'T')
            ++conv[fpos];
          else if (qc != 'N')
            ++err[fpos];
        }
      }
    }
    else {
      if (cigar_eats_query(op)) qpos += n;
      if (cigar_eats_ref(op)) rpos += n;
      if (cigar_eats_frag(op)) fpos += n;
    }
  }

  assert(qpos == get_l_qseq(aln));
}

static void
count_states_neg(const bool INCLUDE_CPGS, const string &chrom,
                 const bam_rec &aln, vector<size_t> &unconv,
                 vector<size_t> &conv, vector<size_t> &err, size_t &hanging) {
  /* iterate backward over query/read positions but forward over
     reference and fragment positions */
  const auto seq = bam_get_seq(aln);
  const auto beg_cig = bam_get_cigar(aln);
  const auto end_cig = beg_cig + get_n_cigar(aln);
  auto rpos = get_pos(aln);
  auto qpos = get_l_qseq(aln);
  auto fpos = 0;

  const decltype(rpos) chrom_lim = chrom.size() - 1;

  for (auto c_itr = beg_cig; c_itr != end_cig; ++c_itr) {
    const auto op = bam_cigar_op(*c_itr);
    const auto n = bam_cigar_oplen(*c_itr);
    if (cigar_eats_ref(op) && cigar_eats_query(op)) {
      const decltype(qpos) end_qpos = qpos - n;
      for (; qpos > end_qpos; --qpos, ++rpos, ++fpos) {
        if (rpos > chrom_lim) ++hanging;
        if (is_guanine(chrom[rpos]) &&
            (rpos == 0 || !is_cytosine(chrom[rpos - 1]) || INCLUDE_CPGS)) {
          const auto qc = seq_nt16_str[bam_seqi(seq, qpos - 1)];
          if (qc == 'C')
            ++unconv[fpos];
          else if (qc == 'T')
            ++conv[fpos];
          else if (qc != 'N')
            ++err[fpos];
        }
      }
    }
    else {
      if (cigar_eats_query(op)) qpos -= n;
      if (cigar_eats_ref(op)) rpos += n;
      if (cigar_eats_frag(op)) fpos += n;
    }
  }
  assert(qpos == 0);
}

static void
write_output(const string &outfile, const vector<size_t> &ucvt_count_p,
             const vector<size_t> &cvt_count_p,
             const vector<size_t> &ucvt_count_n,
             const vector<size_t> &cvt_count_n, const vector<size_t> &err_p,
             const vector<size_t> &err_n) {
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
  if (!out) throw dnmt_error("failed to open output file");

  out << "OVERALL CONVERSION RATE = "
      << static_cast<double>(total_cvt) / (total_cvt + total_ucvt) << endl
      << "POS CONVERSION RATE = "
      << static_cast<double>(pos_cvt) / (pos_cvt + pos_ucvt) << '\t'
      << std::fixed << static_cast<size_t>(pos_cvt + pos_ucvt) << endl
      << "NEG CONVERSION RATE = "
      << static_cast<double>(neg_cvt) / (neg_cvt + neg_ucvt) << '\t'
      << std::fixed << static_cast<size_t>(neg_cvt + neg_ucvt) << endl;

  // clang-format off
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
  // clang-format on

  // Figure out how many positions to print in the output, capped at 1000
  size_t output_len = (ucvt_count_p.size() > 1000) ? 1000 : ucvt_count_p.size();

  while (output_len > 0 &&
         (ucvt_count_p[output_len - 1] + cvt_count_p[output_len - 1] +
            ucvt_count_n[output_len - 1] + cvt_count_n[output_len - 1] ==
          0))
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
        << static_cast<double>(cvt_count_p[i]) / max(size_t(1ul), total_p)
        << '\t';

    out.precision(precision_val);
    out << total_n << '\t' << cvt_count_n[i] << '\t'
        << static_cast<double>(cvt_count_n[i]) / max(size_t(1ul), total_n)
        << '\t';

    const double total_cvt = cvt_count_p[i] + cvt_count_n[i];
    out.precision(precision_val);
    out << static_cast<size_t>(total_valid) << '\t'
        << cvt_count_p[i] + cvt_count_n[i] << '\t'
        << total_cvt / max(1ul, total_valid) << '\t';

    const double total_err = err_p[i] + err_n[i];
    out.precision(precision_val);
    const size_t total = total_valid + err_p[i] + err_n[i];
    out << err_p[i] + err_n[i] << '\t' << static_cast<size_t>(total) << '\t'
        << total_err / max(1ul, total) << endl;
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
    size_t n_threads = 1;

    string chroms_file;
    string outfile;
    string seq_to_use;  // use only this chrom/sequence in the analysis

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]),
                           "Program to compute the "
                           "BS conversion rate from BS-seq "
                           "reads mapped to a genome",
                           "-c <chroms> <mapped-reads>");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)",
                      false, outfile);
    opt_parse.add_opt("chrom", 'c', "File of chromosome sequences (FASTA)",
                      true, chroms_file);
    opt_parse.add_opt("all", 'N', "count all Cs (including CpGs)", false,
                      INCLUDE_CPGS);
    opt_parse.add_opt("seq", '\0', "use only this sequence (e.g. chrM)", false,
                      seq_to_use);
    opt_parse.add_opt("a-rich", 'A', "reads are A-rich", false,
                      reads_are_a_rich);
    opt_parse.add_opt("threads", 't', "number of threads", false, n_threads);
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
    const string bam_file = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/

    vector<string> chroms;
    vector<string> names;
    read_fasta_file_short_names(chroms_file, names, chroms);
    for (auto &&chrom : chroms)
      std::for_each(begin(chrom), end(chrom),
                    [](const char c) { return std::toupper(c); });

    if (VERBOSE)
      cerr << "[n chroms in reference: " << chroms.size() << "]" << endl;

    bamxx::bam_tpool tp(n_threads);

    bamxx::bam_in hts(bam_file);
    if (!hts) throw dnmt_error("failed to open input file: " + bam_file);
    bamxx::bam_header hdr(hts);
    if (!hdr) throw dnmt_error("failed to read header");

    if (n_threads > 1)
      tp.set_io(hts);

    // map the bam header index for each "target" to a sequence in the
    // reference genome
    unordered_map<int32_t, size_t> chrom_lookup;
    size_t chrom_idx_to_use = numeric_limits<size_t>::max();
    for (size_t i = 0; i < chroms.size(); ++i) {
      if (names[i] == seq_to_use) chrom_idx_to_use = i;
      chrom_lookup.insert({sam_hdr_name2tid(hdr.h, names[i].data()), i});
    }

    vector<size_t> unconv_count_pos(output_size, 0ul);
    vector<size_t> conv_count_pos(output_size, 0ul);
    vector<size_t> unconv_count_neg(output_size, 0ul);
    vector<size_t> conv_count_neg(output_size, 0ul);
    vector<size_t> err_pos(output_size, 0ul);
    vector<size_t> err_neg(output_size, 0ul);

    int32_t current_tid = -1;
    size_t chrom_idx = numeric_limits<size_t>::max();
    size_t hanging = 0;

    bool use_this_chrom = seq_to_use.empty();

    bam_rec aln;
    unordered_set<int32_t> chroms_seen;

    while (hts.read(hdr, aln)) {

      if (reads_are_a_rich) flip_conversion(aln);

      // get the correct chrom if it has changed
      if (get_tid(aln) != current_tid) {
        const int32_t the_tid = get_tid(aln);

        // make sure all reads from same chrom are contiguous in the file
        if (chroms_seen.find(the_tid) != end(chroms_seen))
          throw runtime_error("chroms out of order in mapped reads file");

        current_tid = the_tid;

        chroms_seen.insert(the_tid);

        auto chrom_itr = chrom_lookup.find(the_tid);
        if (chrom_itr == end(chrom_lookup))
          throw runtime_error("could not find chrom: " + the_tid);

        chrom_idx = chrom_itr->second;

        if (VERBOSE) cerr << "processing " << names[chrom_idx] << endl;

        use_this_chrom = seq_to_use.empty() || chrom_idx == chrom_idx_to_use;
      }

      if (use_this_chrom) {
        // do the work for this mapped read
        if (bam_is_rev(aln))
          count_states_neg(INCLUDE_CPGS, chroms[chrom_idx], aln,
                           unconv_count_neg, conv_count_neg, err_neg, hanging);
        else
          count_states_pos(INCLUDE_CPGS, chroms[chrom_idx], aln,
                           unconv_count_pos, conv_count_pos, err_pos, hanging);
      }
    }
    write_output(outfile, unconv_count_pos, conv_count_pos, unconv_count_neg,
                 conv_count_neg, err_pos, err_neg);

    if (hanging > 0)  // some overhanging reads
      cerr << "Warning: hanging reads detected at chrom ends "
           << "(N=" << hanging << ")" << endl
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
