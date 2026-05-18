/* bsrate: a program for determining the rate of bisulfite conversion in a
 * bisulfite sequencing experiment
 *
 * Copyright (C) 2014-2026 University of Southern California and
 *                         Andrew D. Smith
 *
 * Authors: Andrew D. Smith and Guilherme Sena
 *
 * This program is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 */

#include "OptionParser.hpp"
#include "bam_record_utils.hpp"
#include "bsutils.hpp"
#include "dnmt_error.hpp"
#include "smithlab_os.hpp"

#include <bamxx.hpp>

#include <htslib/sam.h>

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

static const auto as_frac = [](const auto a, const auto b) {
  return static_cast<double>(a) /
         static_cast<double>(std::max(b, static_cast<decltype(b)>(1)));
};

// NOLINTBEGIN(*-narrowing-conversions,*-pointer-arithmetic)

struct bsrate_summary {
  // converted_count_trich is the number of nucleotides covering a
  // cytosine in the reference that show a thymine in the read, and
  // for reads mapping to the trich strand.
  std::uint64_t converted_count_trich{};

  // total_count_trich is the number of nucleotides covering a
  // cytosine in the reference that show either a cytosine or a
  // thymine in the read, and for reads mapping to the trich
  // strand.
  std::uint64_t total_count_trich{};

  // converted_count_arich is the number of nucleotides covering a
  // cytosine in the reference that show a thymine in the read, and
  // for reads mapping to the arich strand.
  std::uint64_t converted_count_arich{};

  // total_count_arich is the number of nucleotides covering a
  // cytosine in the reference that show either a cytosine or a
  // thymine in the read, and for reads mapping to the arich
  // strand.
  std::uint64_t total_count_arich{};

  // error_count_trich is the number of nucleotides covering a cytosine in
  // the reference shows either an A or a G in a read mapping on the trich
  // strand.
  std::uint64_t error_count_trich{};

  // error_count_arich is the number of nucleotides covering a cytosine in
  // the reference shows either an A or a G in a read mapping on the arich
  // strand.
  std::uint64_t error_count_arich{};

  // total_count is equal to the sum of total_count_trich and
  // total_count_arich
  [[nodiscard]] std::uint64_t
  total_count() const {
    return total_count_trich + total_count_arich;
  }

  // converted_count is equal to the sum of converted_count_trich
  // and converted_count_arich.
  [[nodiscard]] std::uint64_t
  converted_count() const {
    return converted_count_trich + converted_count_arich;
  }

  // bisulfite_conversion_rate_trich is equal to converted_count_trich
  // divided by total_count_trich, a value that is always between 0 and
  // 1. When total_count_trich is 0, then
  // bisulfite_conversion_rate_trich is given a value of 0.
  [[nodiscard]] double
  bisulfite_conversion_rate_trich() const {
    return as_frac(converted_count_trich, total_count_trich);
  }

  // bisulfite_conversion_rate_arich is equal to converted_count_arich
  // divided by total_count_arich, a value that is always between 0 and
  // 1. When total_count_arich is 0, then
  // bisulfite_conversion_rate_arich is given a value of 0.
  [[nodiscard]] double
  bisulfite_conversion_rate_arich() const {
    return as_frac(converted_count_arich, total_count_arich);
  }

  // bisulfite_conversion_rate is equal to converted_count divided by
  // total_count, a value that is always between 0 and 1. When
  // total_count is 0, then bisulfite_conversion_rate is given a value
  // of 0.
  [[nodiscard]] double
  bisulfite_conversion_rate() const {
    return as_frac(converted_count(), total_count());
  }

  // error_count is equal to the sum of error_count_trich and
  // error_count_arich
  [[nodiscard]] std::uint64_t
  error_count() const {
    return error_count_trich + error_count_arich;
  }

  // valid_count_trich is the number of nucleotides in T-rich reads covering a
  // cytosine in the reference with any nucleotide that is not an N in the read.
  [[nodiscard]] std::uint64_t
  valid_count_trich() const {
    return total_count_trich + error_count_trich;
  }

  // valid_count_arich is the number of nucleotides in A-rich reads covering a
  // guanine in the reference with any nucleotide that is not an N in the read.
  [[nodiscard]] std::uint64_t
  valid_count_arich() const {
    return total_count_arich + error_count_arich;
  }

  // valid_count is the number of nucleotides covering a cytosine in the
  // reference with T-rich reads or a guanine with A-rich reads and showing any
  // nucleotide that is not an N in the read.
  [[nodiscard]] std::uint64_t
  valid_count() const {
    return total_count() + error_count();
  }

  // error_rate_trich is equal to error_count_trich divided by
  // valid_count_trich, and is a value that is always between 0 and 1. When
  // valid_count_trich is 0, then error_rate is given a value of 0.
  [[nodiscard]] double
  error_rate_trich() const {
    return as_frac(error_count_trich, valid_count_trich());
  }

  // error_rate_arich is equal to error_count_arich divided by
  // valid_count_arich, and is a value that is always between 0 and 1. When
  // valid_count_arich is 0, then error_rate is given a value of 0.
  [[nodiscard]] double
  error_rate_arich() const {
    return as_frac(error_count_arich, valid_count_arich());
  }

  // error_rate is equal to error_count divided by valid_count, and is
  // a value that is always between 0 and 1. When valid_count is 0,
  // then error_rate is given a value of 0.
  [[nodiscard]] double
  error_rate() const {
    return as_frac(error_count(), valid_count());
  }

  void
  update_trich(const char nt) {
    if (nt == 'C' || nt == 'T') {
      ++total_count_trich;
      converted_count_trich += (nt == 'T');
    }
    else if (nt != 'N')
      ++error_count_trich;
  }

  void
  update_arich(const char nt) {
    if (nt == 'G' || nt == 'A') {
      ++total_count_arich;
      converted_count_arich += (nt == 'A');
    }
    else if (nt != 'N')
      ++error_count_arich;
  }

  bsrate_summary &
  operator+=(const bsrate_summary &rhs) {
    converted_count_trich += rhs.converted_count_trich;
    total_count_trich += rhs.total_count_trich;
    converted_count_arich += rhs.converted_count_arich;
    total_count_arich += rhs.total_count_arich;
    error_count_trich += rhs.error_count_trich;
    error_count_arich += rhs.error_count_arich;
    return *this;
  }

  std::string
  tostring_as_row() const {
    static constexpr auto precision_val = 5u;
    std::ostringstream oss;
    oss.precision(precision_val);
    oss.setf(std::ios_base::fixed, std::ios_base::floatfield);
    // clang-format off
    oss << total_count_trich << '\t'
        << converted_count_trich << '\t'
        << bisulfite_conversion_rate_trich() << '\t';
    oss << total_count_arich << '\t'
        << converted_count_arich << '\t'
        << bisulfite_conversion_rate_arich() << '\t';
    oss << total_count() << '\t'
        << converted_count() << '\t'
        << bisulfite_conversion_rate() << '\t';
    oss << error_count_trich << '\t'
        << valid_count_trich() << '\t'
        << error_rate_trich() << '\t';
    oss << error_count_arich << '\t'
        << valid_count_arich() << '\t'
        << error_rate_arich() << '\t';
    oss << error_count() << '\t'
        << valid_count() << '\t'
        << error_rate();
    // clang-format on
    return oss.str();
  }

  std::string
  tostring_as_yaml_list(const std::uint32_t position) const {
    static constexpr auto precision_val = 5u;
    std::ostringstream oss;
    oss.precision(precision_val);
    oss.setf(std::ios_base::fixed, std::ios_base::floatfield);
    // clang-format off
    oss << "  - base: " << position << '\n'
        << "    ptot: " << total_count_trich << '\n'
        << "    pconv: " << converted_count_trich << '\n'
        << "    prate: " << bisulfite_conversion_rate_trich() << '\n'
        << "    ntot: " << total_count_arich << '\n'
        << "    nconv: " << converted_count_arich << '\n'
        << "    nrate: " << bisulfite_conversion_rate_arich() << '\n'
        << "    bthtot: " << total_count() << '\n'
        << "    bthconv: " << converted_count() << '\n'
        << "    bthrate: " << bisulfite_conversion_rate() << '\n'
        << "    err_trich: " << error_count_trich << '\n'
        << "    all_trich: " << valid_count_trich() << '\n'
        << "    errrate_trich: " << error_rate_trich() << '\n'
        << "    err_arich: " << error_count_arich << '\n'
        << "    all_arich: " << valid_count_arich() << '\n'
        << "    errrate_arich: " << error_rate_arich() << '\n'
        << "    err: " << error_count() << '\n'
        << "    all: " << valid_count() << '\n'
        << "    errrate: " << error_rate() << '\n';
    // clang-format on
    return oss.str();
  }

  std::string
  tostring() const {
    static constexpr auto sep = ": ";
    std::ostringstream oss;
    oss << "converted_count_trich" << sep << converted_count_trich << '\n';
    oss << "total_count_trich" << sep << total_count_trich << '\n';
    oss << "bisulfite_conversion_rate_trich" << sep
        << bisulfite_conversion_rate_trich() << '\n';
    oss << "converted_count_arich" << sep << converted_count_arich << '\n';
    oss << "total_count_arich" << sep << total_count_arich << '\n';
    oss << "bisulfite_conversion_rate_arich" << sep
        << bisulfite_conversion_rate_arich() << '\n';
    oss << "converted_count" << sep << converted_count() << '\n';
    oss << "total_count" << sep << total_count() << '\n';
    oss << "bisulfite_conversion_rate" << sep << bisulfite_conversion_rate()
        << '\n';

    oss << "error_count_trich" << sep << error_count_trich << '\n';
    oss << "valid_count_trich" << sep << valid_count_trich() << '\n';
    oss << "error_rate_trich" << sep << error_rate_trich() << '\n';

    oss << "error_count_arich" << sep << error_count_arich << '\n';
    oss << "valid_count_arich" << sep << valid_count_arich() << '\n';
    oss << "error_rate_arich" << sep << error_rate_arich() << '\n';

    oss << "error_count" << sep << error_count() << '\n';
    oss << "valid_count" << sep << valid_count() << '\n';
    oss << "error_rate" << sep << error_rate();

    return oss.str();
  }
};

inline bsrate_summary
operator+(bsrate_summary lhs, const bsrate_summary &rhs) {
  lhs += rhs;
  return lhs;
}

static std::pair<std::uint32_t, std::uint32_t>
count_states_trich(const bool INCLUDE_CPGS, const std::string &chrom,
                   const bamxx::bam_rec &aln,
                   std::vector<bsrate_summary> &summaries,
                   std::size_t &hanging) {
  std::uint32_t n_conv{};
  std::uint32_t n_uconv{};

  // iterate through reference, query/read and fragment
  const auto seq = bam_get_seq(aln);
  const auto beg_cig = bam_get_cigar(aln);
  const auto end_cig = beg_cig + get_n_cigar(aln);
  auto rpos = get_pos(aln);
  auto qpos = 0;
  auto fpos = 0;

  const auto is_cpg = [&chrom](const std::uint64_t p) {
    return p < std::size(chrom) && chrom[p + 1] == 'G';
  };

  const decltype(rpos) chrom_lim =
    std::size(chrom) > 0 ? std::size(chrom) - 1 : 0;

  for (auto c_itr = beg_cig; c_itr != end_cig; ++c_itr) {
    const auto op = bam_cigar_op(*c_itr);
    const auto n = bam_cigar_oplen(*c_itr);
    if (cigar_eats_ref(op) && cigar_eats_query(op)) {
      const decltype(qpos) end_qpos = qpos + n;
      for (; qpos < end_qpos; ++qpos, ++rpos, ++fpos) {
        if (rpos > chrom_lim)
          ++hanging;
        if (chrom[rpos] == 'C' && (!is_cpg(rpos) || INCLUDE_CPGS)) {
          const auto nt = seq_nt16_str[bam_seqi(seq, qpos)];
          summaries[fpos].update_trich(nt);
          n_conv += (nt == 'T');
          n_uconv += (nt == 'C');
        }
      }
    }
    else {
      if (cigar_eats_query(op))
        qpos += n;
      if (cigar_eats_ref(op))
        rpos += n;
      if (cigar_eats_frag(op))
        fpos += n;
    }
  }
  assert(qpos == get_l_qseq(aln));
  return {n_conv, n_conv + n_uconv};
}

static std::pair<std::uint32_t, std::uint32_t>
count_states_arich(const bool INCLUDE_CPGS, const std::string &chrom,
                   const bamxx::bam_rec &aln,
                   std::vector<bsrate_summary> &summaries,
                   std::size_t &hanging) {
  std::uint32_t n_conv{};
  std::uint32_t n_uconv{};

  // iterate through reference, query/read and fragment
  const auto seq = bam_get_seq(aln);
  const auto beg_cig = bam_get_cigar(aln);
  const auto end_cig = beg_cig + get_n_cigar(aln);
  auto rpos = get_pos(aln);
  auto qpos = 0;
  auto fpos = 0;

  const auto is_cpg = [&chrom](const std::uint64_t p) {
    return p > 0ul && chrom[p - 1] == 'C';
  };

  const decltype(rpos) chrom_lim =
    std::size(chrom) > 0 ? std::size(chrom) - 1 : 0;

  for (auto c_itr = beg_cig; c_itr != end_cig; ++c_itr) {
    const auto op = bam_cigar_op(*c_itr);
    const auto n = bam_cigar_oplen(*c_itr);
    if (cigar_eats_ref(op) && cigar_eats_query(op)) {
      const decltype(qpos) end_qpos = qpos + n;
      for (; qpos < end_qpos; ++qpos, ++rpos, ++fpos) {
        if (rpos > chrom_lim)
          ++hanging;
        if (chrom[rpos] == 'G' && (!is_cpg(rpos) || INCLUDE_CPGS)) {
          const auto nt = seq_nt16_str[bam_seqi(seq, qpos)];
          summaries[fpos].update_arich(nt);
          n_conv += (nt == 'A');
          n_uconv += (nt == 'G');
        }
      }
    }
    else {
      if (cigar_eats_query(op))
        qpos += n;
      if (cigar_eats_ref(op))
        rpos += n;
      if (cigar_eats_frag(op))
        fpos += n;
    }
  }
  assert(qpos == get_l_qseq(aln));
  return {n_conv, n_conv + n_uconv};
}

static void
write_output(const std::string &outfile,
             const std::vector<bsrate_summary> &summaries) {
  static constexpr auto max_output_len = 1000ul;  // cap on frag len to report
  std::ofstream of;
  if (!outfile.empty())
    of.open(outfile.c_str());
  std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
  if (!out)
    throw dnmt_error("failed to open output file");

  const bsrate_summary overall_summary =
    std::reduce(std::cbegin(summaries), std::cend(summaries));

  out << "OVERALL CONVERSION RATE = "
      << overall_summary.bisulfite_conversion_rate() << '\n'
      << "T RICH CONVERSION RATE = "
      << overall_summary.bisulfite_conversion_rate_trich() << '\t'
      << overall_summary.total_count_trich << '\n'
      << "A RICH CONVERSION RATE = "
      << overall_summary.bisulfite_conversion_rate_arich() << '\t'
      << overall_summary.total_count_arich << '\n';

  // clang-format off
  out << "BASE" << '\t'
      << "T_TOT" << '\t'
      << "T_CONV" << '\t'
      << "T_RATE" << '\t'
      << "A_TOT" << '\t'
      << "A_CONV" << '\t'
      << "A_RATE" << '\t'
      << "BTHTOT" << '\t'
      << "BTHCONV" << '\t'
      << "BTHRATE" << '\t'
      << "ERR_T" << '\t'
      << "ALL_T" << '\t'
      << "ERATE_T"  << '\t'
      << "ERR_A" << '\t'
      << "ALL_A" << '\t'
      << "ERATE_A"  << '\t'
      << "ERR" << '\t'
      << "ALL" << '\t'
      << "ERATE"  << '\n';
  // clang-format on

  // figure out how many positions to print in the output
  auto output_len = std::min(std::size(summaries), max_output_len);
  while (output_len > 0 && summaries[output_len - 1].total_count() == 0)
    --output_len;

  for (auto i = 0u; i < output_len; ++i)
    out << (i + 1) << '\t' << summaries[i].tostring_as_row() << '\n';
}

static void
write_output_yaml(const std::string &outfile,
                  const std::vector<bsrate_summary> &summaries) {
  static constexpr auto max_output_len = 1000ul;  // cap on frag len to report
  std::ofstream of;
  if (!outfile.empty())
    of.open(outfile.c_str());
  std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
  if (!out)
    throw dnmt_error("failed to open output file");

  const bsrate_summary overall_summary =
    std::reduce(std::cbegin(summaries), std::cend(summaries));
  // overall_summary.set_values();
  out << "overall_conversion_rate: "
      << overall_summary.bisulfite_conversion_rate() << '\n'
      << "trich_conversion_rate: "
      << overall_summary.bisulfite_conversion_rate_trich() << '\n'
      << "trich_count: " << overall_summary.total_count_trich << '\n'
      << "arich_conversion_rate: "
      << overall_summary.bisulfite_conversion_rate_arich() << '\n'
      << "arich_count: " << overall_summary.total_count_arich << '\n';

  // figure out how many positions to print in the output
  auto output_len = std::min(std::size(summaries), max_output_len);
  while (output_len > 0 && summaries[output_len - 1].total_count() == 0)
    --output_len;

  out << "rows:\n";
  for (auto i = 0u; i < output_len; ++i)
    out << summaries[i].tostring_as_yaml_list(i + 1);
}

static void
write_summary(const std::string &summary_file,
              std::vector<bsrate_summary> &summaries) {
  const auto s = std::reduce(std::cbegin(summaries), std::cend(summaries));
  std::ofstream out(summary_file);
  if (!out)
    throw dnmt_error("failed to open file: " + summary_file);
  out << s.tostring() << '\n';
}

template <typename T>
static inline void
update_per_read_stats(const std::pair<T, T> &x,
                      std::vector<std::vector<T>> &tab) {
  if (x.second < std::size(tab))
    ++tab[x.second][x.first];
}

static inline std::vector<double>
format_histogram(const std::vector<std::vector<std::uint32_t>> &tab,
                 const std::size_t n_hist_bins) {
  static constexpr auto epsilon = 1e-6;
  std::vector<double> hist(n_hist_bins, 0.0);
  for (std::size_t i = 1; i < std::size(tab); ++i) {
    const double denom = i + epsilon;
    for (std::size_t j = 0; j <= i; ++j) {
      const double frac = j / denom;
      const auto bin_id = std::floor(frac * n_hist_bins);
      assert(bin_id < std::size(hist));
      hist[bin_id] += tab[i][j];
    }
  }
  return hist;
}

static void
write_hanging_read_message(std::ostream &out, const std::size_t n_hanging) {
  out << "Warning: hanging reads detected at chromosome ends "
      << "(N=" << n_hanging << ")." << '\n'
      << "High numbers of hanging reads suggest inconsistent " << '\n'
      << "reference genomes between stages of analysis." << '\n'
      << "This is likely to result in analysis errors." << '\n';
}

template <typename T>
static void
write_per_read_histogram(const std::vector<std::vector<T>> &tab,
                         const std::size_t n_hist_bins, std::ostream &out) {
  const auto hist = format_histogram(tab, n_hist_bins);
  out << std::fixed;
  // NOLINTNEXTLINE(cert-flp30-c,clang-analyzer-security*)
  for (auto i = 0.0; i < std::size(hist); ++i)
    out << std::setprecision(3) << i / std::size(hist) << '\t'
        << std::setprecision(3) << (i + 1) / std::size(hist) << '\t'
        << std::setprecision(0) << hist[i] << '\n';
}

int
main_bsrate(int argc, char *argv[]) {  // NOLINT(*-avoid-c-arrays)
  try {
    // default number of bins in a histogram when writing per-read conversion
    // histogram
    static constexpr const std::size_t n_hist_bins_default = 20;

    // assumed maximum length of a fragment
    static constexpr const std::size_t output_size = 10000;

    // Assumed maximum cytosines per fragment. Currently the per-read
    // information is collected as counts in a 2D array, so that later
    // features may be able to fit distributions based these counts.
    static constexpr const std::size_t max_cytosine_per_frag = 1000;

    bool VERBOSE = false;
    bool INCLUDE_CPGS = false;
    bool output_yaml = false;
    bool report_per_read = false;
    std::size_t n_threads = 1;
    std::size_t n_hist_bins = n_hist_bins_default;

    std::string chroms_file;
    std::string summary_file;
    std::string outfile;
    std::string seq_to_use;  // use only this chrom/sequence in the analysis

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(argv[0],  // NOLINT(*-pointer-arithmetic)
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
    opt_parse.add_opt("per-read", 'p', "report per-read conversion to terminal",
                      false, report_per_read);
    opt_parse.add_opt("bins", 'b', "number of bins for per-read", false,
                      n_hist_bins);
    opt_parse.add_opt("summary", 'S', "summary file name", false, summary_file);
    opt_parse.add_opt("yaml", 'y', "output in yaml format", false, output_yaml);
    opt_parse.add_opt("threads", 't', "number of threads", false, n_threads);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    std::vector<std::string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      std::cerr << opt_parse.help_message() << '\n'
                << opt_parse.about_message() << '\n';
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      std::cerr << opt_parse.about_message() << '\n';
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      std::cerr << opt_parse.option_missing_message() << '\n';
      return EXIT_SUCCESS;
    }
    if (std::size(leftover_args) != 1) {
      std::cerr << opt_parse.help_message() << '\n';
      return EXIT_SUCCESS;
    }
    const std::string bam_file = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/

    std::vector<std::string> chroms;
    std::vector<std::string> names;
    read_fasta_file_short_names(chroms_file, names, chroms);
    for (auto &&chrom : chroms)
      std::for_each(std::begin(chrom), std::end(chrom),
                    [](const char c) { return std::toupper(c); });

    if (VERBOSE)
      std::cerr << "[n chroms in reference: " << std::size(chroms) << "]"
                << '\n';

    bamxx::bam_tpool tp(n_threads);

    bamxx::bam_in hts(bam_file);
    if (!hts)
      throw dnmt_error("failed to open input file: " + bam_file);
    bamxx::bam_header hdr(hts);
    if (!hdr)
      throw dnmt_error("failed to read header");

    if (n_threads > 1)
      tp.set_io(hts);

    // map the bam header index for each "target" to a sequence in the
    // reference genome
    std::unordered_map<std::int32_t, std::size_t> chrom_lookup;
    std::size_t chrom_idx_to_use = std::numeric_limits<std::size_t>::max();
    for (auto i = 0u; i < std::size(chroms); ++i) {
      if (names[i] == seq_to_use)
        chrom_idx_to_use = i;
      chrom_lookup.emplace(sam_hdr_name2tid(hdr.h, names[i].data()), i);
    }

    std::vector<std::vector<std::uint32_t>> per_read_counts(
      max_cytosine_per_frag,
      std::vector<std::uint32_t>(max_cytosine_per_frag, 0));

    std::vector<bsrate_summary> summaries(output_size);

    std::int32_t current_tid = -1;
    std::size_t chrom_idx = std::numeric_limits<std::size_t>::max();
    std::size_t hanging = 0;

    bool use_this_chrom = seq_to_use.empty();

    bamxx::bam_rec aln;
    std::unordered_set<std::int32_t> chroms_seen;

    while (hts.read(hdr, aln)) {
      const std::int32_t the_tid = get_tid(aln);
      if (the_tid == -1)
        continue;

      const auto original_is_arich = aln_is_a_rich(aln);

      // get the correct chrom if it has changed
      if (the_tid != current_tid) {
        // make sure all reads from same chrom are contiguous in the file
        if (chroms_seen.find(the_tid) != end(chroms_seen))
          throw std::runtime_error("chroms out of order in mapped reads file");

        current_tid = the_tid;

        chroms_seen.insert(the_tid);

        auto chrom_itr = chrom_lookup.find(the_tid);
        if (chrom_itr == end(chrom_lookup))
          throw std::runtime_error("could not find chrom: " +
                                   std::to_string(the_tid));

        chrom_idx = chrom_itr->second;

        if (VERBOSE)
          std::cerr << "processing " << names[chrom_idx] << '\n';

        use_this_chrom = seq_to_use.empty() || chrom_idx == chrom_idx_to_use;
      }

      if (use_this_chrom) {
        // do the work for the current mapped read
        const auto conv_result =
          original_is_arich
            ? count_states_arich(INCLUDE_CPGS, chroms[chrom_idx], aln,
                                 summaries, hanging)
            : count_states_trich(INCLUDE_CPGS, chroms[chrom_idx], aln,
                                 summaries, hanging);
        if (report_per_read)
          update_per_read_stats(conv_result, per_read_counts);
      }
    }

    if (output_yaml)
      write_output_yaml(outfile, summaries);
    else
      write_output(outfile, summaries);

    if (hanging > 0)
      write_hanging_read_message(std::cerr, hanging);

    if (!summary_file.empty())
      write_summary(summary_file, summaries);

    if (report_per_read)
      write_per_read_histogram(per_read_counts, n_hist_bins, std::cout);
  }
  catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

// NOLINTEND(*-narrowing-conversions,*-pointer-arithmetic)
