/* guessprotocol: a program for guessing whether a wgbs protocol is
 * wgbs, pbat or random pbat
 *
 * Copyright (C) 2019-2023 Andrew D. Smith
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

#include "OptionParser.hpp"
#include "numerical_utils.hpp"

#include <bamxx.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// NOLINTBEGIN(*-pointer-arithmetic,*-constant-array-index,*-narrowing-conversions,cert-err09-cpp,cert-err61-cpp)

// clang-format off
constexpr std::array<int, 256> nuc_to_idx = {
 /*  0*/  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 /* 16*/  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 /* 32*/  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 /* 48*/  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 /* 64*/  4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
 /* 80*/  4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 /* 96*/  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 /*112*/  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 /*128*/  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
};
// clang-format on

struct nucleotide_model {
  std::vector<double> pr{};
  std::vector<double> lpr{};
  double bisulfite_conversion_rate{};
  bool is_t_rich{};

  nucleotide_model(const std::vector<double> &bc, const double conv_rate,
                   const bool itr) :
    pr{bc}, bisulfite_conversion_rate{conv_rate}, is_t_rich{itr} {
    auto nuc_from = is_t_rich ? 1 : 2;
    auto nuc_to = is_t_rich ? 3 : 0;
    pr[nuc_to] += bisulfite_conversion_rate * pr[nuc_from];
    pr[nuc_from] *= (1.0 - bisulfite_conversion_rate);
    assert(std::reduce(std::cbegin(pr), std::cend(pr), 0.0) == 1.0);
    lpr.resize(std::size(pr));
    transform(std::cbegin(pr), std::cend(pr), std::begin(lpr),
              [](const double x) { return log(x); });
  }

  double
  operator()(const std::string &s) const {
    return std::accumulate(std::cbegin(s), std::cend(s), 0.0,
                           [&](const double x, const char c) {
                             const auto i = nuc_to_idx[static_cast<uint8_t>(c)];
                             return i == 4 ? x : x + lpr[i];
                           });
  };

  std::string
  tostring() const {
    std::ostringstream oss;
    oss << "pr:\n";
    for (auto i : pr)
      oss << i << '\n';
    oss << "log pr:\n";
    for (auto i : lpr)
      oss << i << '\n';
    oss << bisulfite_conversion_rate << '\n' << is_t_rich;
    return oss.str();
  }
};

struct guessprotocol_summary {
  static constexpr auto wgbs_cutoff_confident = 0.99;
  static constexpr auto wgbs_cutoff_unconfident = 0.9;
  static constexpr auto rpbat_cutoff_confident_high = 0.8;
  static constexpr auto rpbat_cutoff_confident_low = 0.2;
  static constexpr auto pbat_cutoff_unconfident = 0.1;
  static constexpr auto pbat_cutoff_confident = 0.01;

  // protocol is the guessed protocol (wgbs, pbat or rpbat) based on the
  // content of the reads.
  std::string protocol;
  // confidence indicates the level of confidence in the guess for the
  // protocol.
  std::string confidence;
  // layout indicates whether the reads are paired or single-ended.
  std::string layout;
  // n_reads_wgbs is the average number of reads (for single-ended reads) or
  // read pairs (for paired reads) where read1 is T-rich.
  double n_reads_wgbs{};
  // n_reads is the number of evaluated reads or read pairs.
  uint64_t n_reads{};
  // wgbs_fraction is the probability that a read (for single-ended reads) or
  // the read1 of a read pair (for paired reads) is T-rich.
  double wgbs_fraction{};

  void
  evaluate() {
    const auto frac = n_reads_wgbs / n_reads;

    // assigning wgbs (near one)
    if (frac > wgbs_cutoff_confident) {
      protocol = "wgbs";
      confidence = "high";
    }
    else if (frac > wgbs_cutoff_unconfident) {
      protocol = "wgbs";
      confidence = "low";
    }
    // assigning pbat (near zero)
    else if (frac < pbat_cutoff_confident) {
      protocol = "pbat";
      confidence = "high";
    }
    else if (frac < pbat_cutoff_unconfident) {
      protocol = "pbat";
      confidence = "low";
    }
    // assigning rpbat (towards middle)
    else if (frac > rpbat_cutoff_confident_low &&
             frac < rpbat_cutoff_confident_high) {
      protocol = "rpbat";
      confidence = "high";
    }
    else {
      protocol = "rpbat";
      confidence = "low";
    }

    wgbs_fraction = frac;
  }

  std::string
  tostring() const {
    std::ostringstream oss;
    oss << "protocol: " << protocol << '\n'
        << "confidence: " << confidence << '\n'
        << "wgbs_fraction: " << wgbs_fraction << '\n'
        << "n_reads_wgbs: " << n_reads_wgbs << '\n'
        << "n_reads: " << n_reads;
    return oss.str();
  }
};

// store each read from one end
struct FASTQRecord {
  std::string name;
  std::string seq;
};

// see if two reads from two ends match to each other (they should
// have the same name)
static bool
mates(const size_t to_ignore_at_end,  // in case names have #0/1 name ends
      const FASTQRecord &a, const FASTQRecord &b) {
  assert(to_ignore_at_end < std::size(a.name));
  return equal(std::cbegin(a.name), std::cend(a.name) - to_ignore_at_end,
               std::cbegin(b.name));
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
  r.name.erase(std::copy_n(std::cbegin(r.name) + 1, nm_sz, std::begin(r.name)),
               std::cend(r.name));

  if (!getline(s, r.seq))
    ec = err_code::bad_seq;

  std::string tmp;
  if (!getline(s, tmp))
    ec = err_code::bad_plus;

  if (!getline(s, tmp))
    ec = err_code::bad_qual;

  if (ec != err_code::none)
    throw error_msg[ec];  // NOLINT(runtime/arrays)

  return s;
}

int
main_guessprotocol(int argc, char *argv[]) {  // NOLINT(*-avoid-c-arrays)
  try {
    static const std::vector<double> human_base_comp = {
      0.295,
      0.205,
      0.205,
      0.295,
    };
    static const std::vector<double> flat_base_comp = {
      0.25,
      0.25,
      0.25,
      0.25,
    };

    constexpr auto description = "guess bisulfite protocol for a library";

    bool verbose = false;
    bool use_human = false;
    std::string outfile;
    size_t reads_to_check = 1000000;  // NOLINT(*-avoid-magic-numbers)
    size_t name_suffix_len = 0;
    double bisulfite_conversion_rate = 0.98;  // NOLINT(*-avoid-magic-numbers)

    namespace fs = std::filesystem;
    const std::string cmd_name = std::filesystem::path(argv[0]).filename();

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(cmd_name, description,
                           "<end1-fastq> [<end2-fastq>]");
    opt_parse.add_opt("nreads", 'n', "number of reads in initial check", false,
                      reads_to_check);
    opt_parse.add_opt("ignore", 'i',
                      "length of read name suffix "
                      "to ignore when matching",
                      false, name_suffix_len);
    opt_parse.add_opt("bisulfite", 'b', "bisulfite conversion rate", false,
                      bisulfite_conversion_rate);
    opt_parse.add_opt("human", 'H', "assume human genome", false, use_human);
    opt_parse.add_opt("output", 'o', "output file name", false, outfile);
    opt_parse.add_opt("verbose", 'v',
                      "report available information during the run", false,
                      verbose);
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
    if (opt_parse.about_requested() || leftover_args.size() > 2) {
      std::cerr << opt_parse.about_message() << '\n';
      return EXIT_SUCCESS;
    }
    const std::vector<std::string> reads_files(leftover_args);
    /****************** END COMMAND LINE OPTIONS *****************/

    auto base_comp = flat_base_comp;
    if (use_human)
      base_comp = human_base_comp;

    nucleotide_model t_rich_model(base_comp, bisulfite_conversion_rate, true);
    nucleotide_model a_rich_model(base_comp, bisulfite_conversion_rate, false);

    guessprotocol_summary summary;
    summary.layout = reads_files.size() == 2 ? "paired" : "single";

    if (verbose) {
      if (reads_files.size() == 2)
        std::cerr << "data layout: "
                  << "paired" << '\n'
                  << "read1 file: " << reads_files.front() << '\n'
                  << "read2 file: " << reads_files.back() << '\n';
      else
        std::cerr << "data layout: "
                  << "single" << '\n'
                  << "read file: " << reads_files.front() << '\n';
      std::cerr << "reads to check: " << reads_to_check << '\n'
                << "read name suffix length: " << name_suffix_len << '\n'
                << "bisulfite conversion: " << bisulfite_conversion_rate
                << '\n';
    }

    if (reads_files.size() == 2) {
      // input: paired-end reads with end1 and end2
      bamxx::bgzf_file in1(reads_files.front(), "r");
      if (!in1)
        throw std::runtime_error("cannot open file: " + reads_files.front());

      bamxx::bgzf_file in2(reads_files.back(), "r");
      if (!in2)
        throw std::runtime_error("cannot open file: " + reads_files.back());

      FASTQRecord r1, r2;
      while (in1 >> r1 && in2 >> r2 && summary.n_reads < reads_to_check) {
        summary.n_reads++;

        if (!mates(name_suffix_len, r1, r2))
          throw std::runtime_error("expected mates: " + r1.name + ", " +
                                   r2.name);

        const double ta = t_rich_model(r1.seq) + a_rich_model(r2.seq);
        const double at = a_rich_model(r1.seq) + t_rich_model(r2.seq);

        const auto prob_read1_t_rich = exp(ta - log_sum_log(ta, at));
        summary.n_reads_wgbs += prob_read1_t_rich;
      }
    }
    else {
      // input: single-end reads
      bamxx::bgzf_file in(reads_files.front(), "r");
      if (!in)
        throw std::runtime_error("cannot open file: " + reads_files.front());

      FASTQRecord r;
      while (in >> r && summary.n_reads < reads_to_check) {
        summary.n_reads++;

        const double t = t_rich_model(r.seq);
        const double a = a_rich_model(r.seq);

        const auto prob_t_rich = exp(t - log_sum_log(t, a));
        summary.n_reads_wgbs += prob_t_rich;
      }
    }

    summary.evaluate();

    if (!outfile.empty()) {
      std::ofstream out(outfile);
      if (!out)
        throw std::runtime_error("failed to open: " + outfile);
      out << summary.tostring() << '\n';
    }
    else
      std::cout << summary.tostring() << '\n';
  }
  catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

// NOLINTEND(*-pointer-arithmetic,*-constant-array-index,*-narrowing-conversions,cert-err09-cpp,cert-err61-cpp)
