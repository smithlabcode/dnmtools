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

#include <bamxx.hpp>

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

using std::accumulate;
using std::cerr;
using std::cout;
using std::endl;
using std::for_each;
using std::max;
using std::numeric_limits;
using std::pair;
using std::runtime_error;
using std::string;
using std::unordered_map;
using std::unordered_set;
using std::vector;

using bamxx::bam_rec;

struct bsrate_summary {
  // converted_count_positive is the number of nucleotides covering a
  // cytosine in the reference that show a thymine in the read, and
  // for reads mapping to the positive strand.
  uint64_t converted_count_positive{};
  // total_count_positive is the number of nucleotides covering a
  // cytosine in the reference that show either a cytosine or a
  // thymine in the read, and for reads mapping to the positive
  // strand.
  uint64_t total_count_positive{};
  // bisulfite_conversion_rate_positive is equal to
  // converted_count_positive divided by total_count_positive, a value
  // that is always between 0 and 1. When total_count_positive is 0,
  // then bisulfite_conversion_rate_positive is given a value of 0.
  double bisulfite_conversion_rate_positive{};

  // converted_count_negative is the number of nucleotides covering a
  // cytosine in the reference that show a thymine in the read, and
  // for reads mapping to the negative strand.
  uint64_t converted_count_negative{};
  // total_count_negative is the number of nucleotides covering a
  // cytosine in the reference that show either a cytosine or a
  // thymine in the read, and for reads mapping to the negative
  // strand.
  uint64_t total_count_negative{};
  // bisulfite_conversion_rate_negative is equal to
  // converted_count_negative divided by total_count_negative, a value
  // that is always between 0 and 1. When total_count_negative is 0,
  // then bisulfite_conversion_rate_negative is given a value of 0.
  double bisulfite_conversion_rate_negative{};

  // converted_count is equal to the sum of converted_count_positive
  // and converted_count_negative.
  uint64_t converted_count{};
  // total_count is equal to the sum of total_count_positive and
  // total_count_negative
  uint64_t total_count{};
  // bisulfite_conversion_rate is equal to converted_count divided by
  // total_count, a value that is always between 0 and 1. When
  // total_count is 0, then bisulfite_conversion_rate is given a value
  // of 0.
  double bisulfite_conversion_rate{};

  // error_count is the number of nucleotides covering a cytosine in
  // the reference shows either an A or a G in the read.
  uint64_t error_count{};
  // valid_count is the number of nucleotides covering a cytosine in
  // the reference shows any nucleotide that is not an N in the read.
  uint64_t valid_count{};
  // error_rate is equal to error_count divided by valid_count, and is
  // a value that is always between 0 and 1. When valid_count is 0,
  // then error_rate is given a value of 0.
  double error_rate{};

  void update_pos(const char nt) {
    if (nt == 'C' || nt == 'T') {
      ++total_count_positive;
      converted_count_positive += (nt == 'T');
    }
    else if (nt != 'N')
      ++error_count;
  }

  void update_neg(const char nt) {
    if (nt == 'C' || nt == 'T') {
      ++total_count_negative;
      converted_count_negative += (nt == 'T');
    }
    else if (nt != 'N')
      ++error_count;
  }

  bsrate_summary &operator+=(const bsrate_summary &rhs) {
    // ADS: the "rates" are set to 0.0 here to ensure that they are
    // computed properly using the full data after any accumulation of
    // integer values has completed.
    converted_count_positive += rhs.converted_count_positive;
    total_count_positive += rhs.total_count_positive;
    bisulfite_conversion_rate_positive = 0.0;

    converted_count_negative += rhs.converted_count_negative;
    total_count_negative += rhs.total_count_negative;
    bisulfite_conversion_rate_negative = 0.0;

    converted_count += rhs.converted_count;
    total_count += rhs.total_count;
    bisulfite_conversion_rate = 0.0;

    error_count += rhs.error_count;
    valid_count += rhs.valid_count;
    error_rate = 0.0;
    return *this;
  }

  void set_values() {
    bisulfite_conversion_rate_positive =
      static_cast<double>(converted_count_positive) /
      max(total_count_positive, 1ul);

    bisulfite_conversion_rate_negative =
      static_cast<double>(converted_count_negative) /
      max(total_count_negative, 1ul);

    converted_count = converted_count_positive + converted_count_negative;
    total_count = total_count_positive + total_count_negative;

    bisulfite_conversion_rate =
      static_cast<double>(converted_count) / max(total_count, 1ul);

    valid_count = total_count + error_count;

    error_rate = static_cast<double>(error_count) / max(valid_count, 1ul);
  }

  string tostring_as_row() const {
    static constexpr auto precision_val = 5u;
    std::ostringstream oss;
    oss.precision(precision_val);
    oss.setf(std::ios_base::fixed, std::ios_base::floatfield);
    // clang-format off
    oss << total_count_positive << '\t'
        << converted_count_positive << '\t'
        << bisulfite_conversion_rate_positive << '\t';
    oss << total_count_negative << '\t'
        << converted_count_negative << '\t'
        << bisulfite_conversion_rate_negative << '\t';
    oss << total_count << '\t'
        << converted_count << '\t'
        << bisulfite_conversion_rate << '\t';
    oss << error_count << '\t'
        << valid_count << '\t'
        << error_rate;
    // clang-format on
    return oss.str();
  }

  string tostring() const {
    static constexpr auto sep = ": ";
    std::ostringstream oss;
    oss << "converted_count_positive" << sep << converted_count_positive
        << '\n';

    oss << "total_count_positive" << sep << total_count_positive << '\n';

    oss << "bisulfite_conversion_rate_positive" << sep
        << bisulfite_conversion_rate_positive << '\n';

    oss << "converted_count_negative" << sep << converted_count_negative
        << '\n';

    oss << "total_count_negative" << sep << total_count_negative << '\n';

    oss << "bisulfite_conversion_rate_negative" << sep
        << bisulfite_conversion_rate_negative << '\n';

    oss << "converted_count" << sep << converted_count << '\n';
    oss << "total_count" << sep << total_count << '\n';

    oss << "bisulfite_conversion_rate" << sep << bisulfite_conversion_rate
        << '\n';

    oss << "error_count" << sep << error_count << '\n';
    oss << "valid_count" << sep << valid_count << '\n';
    oss << "error_rate" << sep << error_rate;

    return oss.str();
  }
};

inline bsrate_summary
operator+(bsrate_summary lhs, const bsrate_summary &rhs) {
  lhs += rhs;
  return lhs;
}

static pair<uint32_t, uint32_t>
count_states_pos(const bool INCLUDE_CPGS, const string &chrom,
                 const bam_rec &aln, vector<bsrate_summary> &summaries,
                 size_t &hanging) {
  uint32_t n_conv = 0, n_uconv = 0;

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
          const auto nt = seq_nt16_str[bam_seqi(seq, qpos)];
          summaries[fpos].update_pos(nt);
          n_conv += (nt == 'T');
          n_uconv += (nt == 'C');
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
  return {n_conv, n_conv + n_uconv};
}

static pair<uint32_t, uint32_t>
count_states_neg(const bool INCLUDE_CPGS, const string &chrom,
                 const bam_rec &aln, vector<bsrate_summary> &summaries,
                 size_t &hanging) {
  uint32_t n_conv = 0, n_uconv = 0;

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
          const auto nt = seq_nt16_str[bam_seqi(seq, qpos - 1)];
          summaries[fpos].update_neg(nt);
          n_conv += (nt == 'T');
          n_uconv += (nt == 'C');
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
  return {n_conv, n_conv + n_uconv};
}

static void
write_output(const string &outfile, const vector<bsrate_summary> &summaries) {
  std::ofstream of;
  if (!outfile.empty()) of.open(outfile.c_str());
  std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
  if (!out) throw dnmt_error("failed to open output file");

  bsrate_summary overall_summary = reduce(cbegin(summaries), cend(summaries));
  overall_summary.set_values();

  out << "OVERALL CONVERSION RATE = "
      << overall_summary.bisulfite_conversion_rate << '\n'
      << "POS CONVERSION RATE = "
      << overall_summary.bisulfite_conversion_rate_positive << '\t'
      << overall_summary.total_count_positive << '\n'
      << "NEG CONVERSION RATE = "
      << overall_summary.bisulfite_conversion_rate_negative << '\t'
      << overall_summary.total_count_negative << '\n';

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

  // figure out how many positions to print in the output, capped at 1000
  auto output_len = std::min(std::size(summaries), 1000ul);
  while (output_len > 0 && summaries[output_len - 1].total_count == 0)
    --output_len;

  for (auto i = 0u; i < output_len; ++i)
    out << (i + 1) << '\t' << summaries[i].tostring_as_row() << '\n';
}

static void
write_summary(const string &summary_file, vector<bsrate_summary> &summaries) {
  auto s = reduce(cbegin(summaries), cend(summaries));
  s.set_values();

  std::ofstream out(summary_file);
  if (!out) throw dnmt_error("failed to open file: " + summary_file);
  out << s.tostring() << '\n';
}

template<typename T> static inline void
update_per_read_stats(const pair<T, T> &x, vector<vector<T>> &tab) {
  if (x.second < std::size(tab)) ++tab[x.second][x.first];
}

static inline vector<double>
format_histogram(const vector<vector<uint32_t>> &tab,
                 const size_t n_hist_bins) {
  static constexpr auto epsilon = 1e-6;
  vector<double> hist(n_hist_bins, 0.0);
  for (size_t i = 1; i < tab.size(); ++i) {
    const double denom = i + epsilon;
    for (size_t j = 0; j <= i; ++j) {
      const double frac = j / denom;
      const auto bin_id = std::floor(frac * n_hist_bins);
      assert(bin_id < hist.size());
      hist[bin_id] += tab[i][j];
    }
  }
  return hist;
}

static void
write_hanging_read_message(std::ostream &out, const size_t n_hanging) {
  out << "Warning: hanging reads detected at chromosome ends "
      << "(N=" << n_hanging << ")." << endl
      << "High numbers of hanging reads suggest inconsistent " << endl
      << "reference genomes between stages of analysis." << endl
      << "This is likely to result in analysis errors." << endl;
}

template<typename T> static void
write_per_read_histogram(const vector<vector<T>> &tab, const size_t n_hist_bins,
                         std::ostream &out) {
  const auto hist = format_histogram(tab, n_hist_bins);
  out << std::fixed;
  for (auto i = 0.0; i < hist.size(); ++i)
    out << std::setprecision(3) << i / hist.size() << '\t'
        << std::setprecision(3) << (i + 1) / hist.size() << '\t'
        << std::setprecision(0) << hist[i] << endl;
}

int
main_bsrate(int argc, const char **argv) {
  try {
    // assumed maximum length of a fragment
    static constexpr const size_t output_size = 10000;

    // Assumed maximum cytosines per fragment. Currently the per-read
    // information is collected as counts in a 2D array, so that later
    // features may be able to fit distributions based these counts.
    static constexpr const size_t max_cytosine_per_frag = 1000;

    bool VERBOSE = false;
    bool INCLUDE_CPGS = false;
    bool reads_are_a_rich = false;
    bool report_per_read = false;
    size_t n_threads = 1;
    size_t n_hist_bins = 20;

    string chroms_file;
    string summary_file;
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
    opt_parse.add_opt("per-read", 'p', "report per-read conversion to terminal",
                      false, report_per_read);
    opt_parse.add_opt("bins", 'b', "number of bins for per-read", false,
                      n_hist_bins);
    opt_parse.add_opt("summary", 'S', "summary file name", false, summary_file);
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

    if (n_threads > 1) tp.set_io(hts);

    // map the bam header index for each "target" to a sequence in the
    // reference genome
    unordered_map<int32_t, size_t> chrom_lookup;
    size_t chrom_idx_to_use = numeric_limits<size_t>::max();
    for (auto i = 0u; i < std::size(chroms); ++i) {
      if (names[i] == seq_to_use) chrom_idx_to_use = i;
      chrom_lookup.emplace(sam_hdr_name2tid(hdr.h, names[i].data()), i);
    }

    vector<vector<uint32_t>> per_read_counts(
      max_cytosine_per_frag, vector<uint32_t>(max_cytosine_per_frag, 0));

    vector<bsrate_summary> summaries(output_size);

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
        // do the work for the current mapped read
        const auto conv_result =
          bam_is_rev(aln) ? count_states_neg(INCLUDE_CPGS, chroms[chrom_idx],
                                             aln, summaries, hanging)
                          : count_states_pos(INCLUDE_CPGS, chroms[chrom_idx],
                                             aln, summaries, hanging);
        if (report_per_read)
          update_per_read_stats(conv_result, per_read_counts);
      }
    }

    // before writing the output we need to ensure the derived
    // quantities in bsrate_summary have been calculated
    for_each(begin(summaries), end(summaries),
             [](bsrate_summary &s) { s.set_values(); });

    write_output(outfile, summaries);

    if (hanging > 0) write_hanging_read_message(cerr, hanging);

    if (!summary_file.empty()) write_summary(summary_file, summaries);

    if (report_per_read)
      write_per_read_histogram(per_read_counts, n_hist_bins, cout);
  }
  catch (const runtime_error &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
