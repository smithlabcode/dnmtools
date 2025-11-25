/* clean-hairpins: a program for identifying and removing hairpin
 * reads from paired-end WGBS/PBAT/RPBAT or RRBS reads.
 *
 * Copyright (C) 2024 Andrew D. Smith
 *
 * Author: Andrew D. Smith
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

#include <bamxx.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// NOLINTBEGIN(*-narrowing-conversions,*-magic-numbers)

// store each read from one end
struct FASTQRecord {
  std::string name;
  std::string seq;
  std::string qual;
};

bamxx::bgzf_file &
operator<<(bamxx::bgzf_file &out, const FASTQRecord &r) {
  // below, the other chars are @, + and four newlines
  static const std::uint32_t other_chars = 6;
  const std::size_t buf_size =
    size(r.name) + size(r.seq) + size(r.qual) + other_chars;
  std::string buffer(buf_size, '\0');
  std::string::iterator b = begin(buffer);
  *b++ = '@';
  b = copy(std::cbegin(r.name), std::cend(r.name), b);
  *b++ = '\n';
  b = copy(std::cbegin(r.seq), std::cend(r.seq), b);
  *b++ = '\n';
  *b++ = '+';
  *b++ = '\n';
  b = copy(std::cbegin(r.qual), std::cend(r.qual), b);
  *b++ = '\n';
  out.write(buffer);
  return out;
}

// Read 4 lines one time from fastq and fill in the FASTQRecord structure
static bamxx::bgzf_file &
operator>>(bamxx::bgzf_file &s, FASTQRecord &r) {
  static constexpr auto n_error_codes = 5u;

  enum err_code : std::uint8_t {
    none,
    bad_name,
    bad_seq,
    bad_plus,
    bad_qual,
  };

  static const std::array<std::runtime_error, n_error_codes> error_msg = {
    std::runtime_error(""),
    std::runtime_error("failed to parse fastq name line"),
    std::runtime_error("failed to parse fastq sequence line"),
    std::runtime_error("failed to parse fastq plus line"),
    std::runtime_error("failed to parse fastq qual line"),
  };

  err_code ec = err_code::none;

  if (!getline(s, r.name))
    return s;

  if (r.name.empty() || r.name[0] != '@')
    ec = err_code::bad_name;

  const auto nm_end = r.name.find_first_of(" \t");
  const auto nm_sz = (nm_end == std::string::npos ? r.name.size() : nm_end) - 1;
  r.name.erase(copy_n(std::cbegin(r.name) + 1, nm_sz, begin(r.name)),
               std::cend(r.name));

  if (!getline(s, r.seq))
    ec = err_code::bad_seq;

  std::string tmp;
  if (!getline(s, tmp))
    ec = err_code::bad_plus;

  if (!getline(s, r.qual))
    ec = err_code::bad_qual;

  if (ec != err_code::none)
    throw error_msg[ec];  // NOLINT

  return s;
}

static inline bool
similar_letters_bisulfite_tc_and_ag(const char a, const char b) {
  return (a != 'N' && a == b) || (a == 'T' && b == 'C') ||
         (a == 'G' && b == 'A');
}

// compare two reads to detect the overlapped region
static double
similarity_both_bisulfite_conversions(const std::string &s1,
                                      const std::string &s2) {
  const std::size_t lim = std::min(size(s1), size(s2));

  std::uint32_t total_letters = 0;
  std::uint32_t matching_letters = 0;
  for (std::size_t i = 0; i < lim; ++i) {
    matching_letters += (similar_letters_bisulfite_tc_and_ag(s1[i], s2[i]));
    total_letters += (valid_base(s1[i]) && valid_base(s2[i]));
  }

  return static_cast<double>(matching_letters) / total_letters;
}

struct hp_summary {
  explicit hp_summary(const double cutoff) : cutoff{cutoff} {}

  // n_reads is the total number of read pairs in the input fastq
  // files.
  std::uint64_t n_reads{};

  // n_good_reads is the total number of read pairs that together have
  // a fraction of good bases above the minimum specified.
  std::uint64_t n_good_reads{};

  // n_hairpin_reads is the number of read pairs identified as being
  // hairpins using the criteria in the "cutoff" variable.
  std::uint64_t n_hairpin_reads{};

  // sum_percent_match_good is the sum of the percent matches between
  // the read ends for reads that do not meet the criteria for hairpin.
  double sum_percent_match_good{};

  // sum_percent_match_good is the sum of the percent matches between
  // the read ends for reads that meet the criteria for hairpin.
  double sum_percent_match_bad{};

  // cutoff is the fraction of matches between the two ends of the
  // read when matching under the assumption that the ends are from a
  // hairpin.
  double cutoff{};

  // mean_percent_match_non_hairpin is the ratio of the
  // sum_percent_match_good over the total non-hairpin reads.
  double mean_percent_match_non_hairpin{};

  // mean_percent_match_hairpin is the ratio of the
  // sum_percent_match_bad over the total hairpin reads.
  double mean_percent_match_hairpin{};

  auto
  assign_values() -> void {
    mean_percent_match_non_hairpin =
      sum_percent_match_good / (n_reads - n_hairpin_reads);
    mean_percent_match_hairpin = sum_percent_match_bad / n_hairpin_reads;
  }

  auto
  tostring() const -> std::string {
    std::ostringstream oss;
    oss << "total_reads_pairs: " << n_reads << '\n'
        << "good_reads_pairs: " << n_good_reads << '\n'
        << "hairpin_read_pairs: " << n_hairpin_reads << '\n'
        << "hairpin_cutoff: " << cutoff << '\n'
        << "sum_percent_match_good: " << sum_percent_match_good << '\n'
        << "mean_percent_match_non_hairpin: " << mean_percent_match_non_hairpin
        << '\n'
        << "sum_percent_match_bad: " << sum_percent_match_bad << '\n'
        << "mean_percent_match_hairpin: " << mean_percent_match_hairpin << '\n';
    return oss.str();
  }
};

static void
write_histogram(const std::string &hist_outfile, std::vector<double> hist) {
  std::ofstream hist_out(hist_outfile);
  if (!hist_out)
    throw std::runtime_error("failed to open file: " + hist_outfile);
  const auto total = std::accumulate(std::cbegin(hist), std::cend(hist), 0.0);
  std::transform(std::cbegin(hist), std::cend(hist), std::begin(hist),
                 [&](const double t) { return t / total; });
  const double increment = 1.0 / size(hist);
  auto idx = 0;
  hist_out << std::fixed;
  hist_out.precision(3);
  // NOLINTBEGIN(clang-analyzer-*,cert-flp30-c)
  for (double offset = 0.0; offset < 1.0; offset += increment)
    hist_out << offset << '\t' << hist[idx++] << '\n';
  // NOLINTEND(clang-analyzer-*,cert-flp30-c)
}

static void
write_statistics(const std::string &filename, hp_summary hps) {
  hps.assign_values();
  std::ofstream out(filename);
  if (!out)
    throw std::runtime_error("failed to open file: " + filename);
  out << hps.tostring();
}

static inline double
fraction_good_bases(const FASTQRecord &a, const FASTQRecord &b) {
  const double a_bases =
    std::count_if(std::cbegin(a.seq), std::cend(a.seq), &valid_base);
  const double b_bases =
    std::count_if(std::cbegin(b.seq), std::cend(b.seq), &valid_base);
  return (a_bases + b_bases) / (size(a.seq) + size(b.seq));
}

struct clean_hairpin {
  double cutoff{0.95};
  // ADS: this was std::uint64_t but g++-14.2.0 on macOS had a problem
  std::size_t n_reads_to_check{std::numeric_limits<std::size_t>::max()};
  double min_good_base_percent{0.5};
  std::uint32_t min_read_length{32};
  std::uint32_t n_hist_bins{20};
  bool invert_output{false};

  hp_summary
  analyze_reads(const std::string &outfile1, const std::string &outfile2,
                bamxx::bgzf_file &in1, bamxx::bgzf_file &in2,
                std::vector<double> &hist) const;
  hp_summary
  analyze_reads(bamxx::bgzf_file &in1, bamxx::bgzf_file &in2,
                std::vector<double> &hist) const;
};

hp_summary
clean_hairpin::analyze_reads(const std::string &outfile1,
                             const std::string &outfile2, bamxx::bgzf_file &in1,
                             bamxx::bgzf_file &in2,
                             std::vector<double> &hist) const {
  // output files for read1 and read2 with hairpins removed
  bamxx::bgzf_file out1(outfile1, "w");
  if (!out1)
    throw std::runtime_error("cannot open output file: " + outfile1);
  bamxx::bgzf_file out2(outfile2, "w");
  if (!out2)
    throw std::runtime_error("cannot open output file: " + outfile2);

  hp_summary hps{cutoff};
  FASTQRecord r1, r2;
  while (hps.n_reads < n_reads_to_check && in1 >> r1 && in2 >> r2) {
    ++hps.n_reads;

    if (fraction_good_bases(r1, r2) < min_good_base_percent ||
        size(r1.seq) < min_read_length || size(r2.seq) < min_read_length)
      continue;

    ++hps.n_good_reads;

    // See if inverted duplicates emerge
    const double percent_match =
      similarity_both_bisulfite_conversions(r1.seq, r2.seq);

    // ADS: need a bitter way to get this bin identifier
    ++hist[std::floor(percent_match * n_hist_bins)];

    if (percent_match > cutoff) {
      hps.sum_percent_match_bad += percent_match;
      hps.n_hairpin_reads++;
      if (invert_output) {
        out1 << r1;
        out2 << r2;
      }
    }
    else {
      hps.sum_percent_match_good += percent_match;
      if (!invert_output) {
        out1 << r1;
        out2 << r2;
      }
    }
  }
  return hps;
}

hp_summary
clean_hairpin::analyze_reads(bamxx::bgzf_file &in1, bamxx::bgzf_file &in2,
                             std::vector<double> &hist) const {
  hp_summary hps{cutoff};
  FASTQRecord r1, r2;
  while (hps.n_reads < n_reads_to_check && in1 >> r1 && in2 >> r2) {
    ++hps.n_reads;

    if (fraction_good_bases(r1, r2) < min_good_base_percent ||
        size(r1.seq) < min_read_length || size(r2.seq) < min_read_length)
      continue;

    ++hps.n_good_reads;

    const double percent_match =
      similarity_both_bisulfite_conversions(r1.seq, r2.seq);
    ++hist[std::floor(percent_match * n_hist_bins)];

    if (percent_match > cutoff) {
      hps.sum_percent_match_bad += percent_match;
      hps.n_hairpin_reads++;
    }
    else
      hps.sum_percent_match_good += percent_match;
  }
  return hps;
}

int
main_clean_hairpins(int argc, char *argv[]) {  // NOLINT(*-avoid-c-arrays)
  static constexpr auto description = "fix and stat invdup/hairpin reads";
  try {
    std::string outfile1;
    std::string outfile2;
    std::string stat_outfile;
    std::string hist_outfile;

    clean_hairpin ch{};

    bool verbose = false;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(argv[0],  // NOLINT(*-pointer-arithmetic)
                           description, "<end1-fastq> <end2-fastq>", true);
    opt_parse.set_show_defaults();
    opt_parse.add_opt("out1", 'o', "output file for read 1", false, outfile1);
    opt_parse.add_opt("out2", 'p', "output file for read 2", false, outfile2);
    opt_parse.add_opt("stat", 's', "stats output file", false, stat_outfile);
    opt_parse.add_opt("nreads", 'n', "number of reads to process", false,
                      ch.n_reads_to_check);
    opt_parse.add_opt("cutoff", 'c', "min percent id to call a hairpin", false,
                      ch.cutoff);
    opt_parse.add_opt("good-bases", 'g', "min percent non-N to consider", false,
                      ch.min_good_base_percent);
    opt_parse.add_opt("length", 'l', "min read length to consider", false,
                      ch.min_read_length);
    opt_parse.add_opt("invert", 'I', "output hairpins instead good pairs",
                      false, ch.invert_output);
    opt_parse.add_opt("hist", 'H', "write a histogram of hairpin matches here",
                      false, hist_outfile);
    opt_parse.add_opt("bins", 'b', "number of histograms bins", false,
                      ch.n_hist_bins);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, verbose);
    std::vector<std::string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      std::cerr << opt_parse.help_message() << '\n'
                << opt_parse.about_message() << '\n';
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      std::cerr << opt_parse.option_missing_message() << '\n';
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested() || leftover_args.size() != 2) {
      std::cerr << opt_parse.about_message() << '\n';
      return EXIT_SUCCESS;
    }
    const std::string reads_file1(leftover_args.front());
    const std::string reads_file2(leftover_args.back());
    /****************** END COMMAND LINE OPTIONS *****************/

    if (outfile1.empty() != outfile2.empty()) {
      // ADS: add message about number of reads to check and output files
      std::cerr << "error: specify both or neither of out1/o and out2/p"
                << '\n';
      return EXIT_FAILURE;
    }

    const bool write_reads = !outfile1.empty();

    std::vector<double> hist(ch.n_hist_bins, 0.0);

    // Input: paired-end reads with end1 and end2
    bamxx::bgzf_file in1(reads_file1, "r");
    if (!in1)
      throw std::runtime_error("cannot open input file: " + reads_file1);
    bamxx::bgzf_file in2(reads_file2, "r");
    if (!in2)
      throw std::runtime_error("cannot open input file: " + reads_file2);

    const hp_summary hps =
      write_reads ? ch.analyze_reads(outfile1, outfile2, in1, in2, hist)
                  : ch.analyze_reads(in1, in2, hist);

    if (!stat_outfile.empty())
      write_statistics(stat_outfile, hps);

    if (!hist_outfile.empty())
      write_histogram(hist_outfile, hist);
  }
  catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

// NOLINTEND(*-narrowing-conversions,*-magic-numbers)
