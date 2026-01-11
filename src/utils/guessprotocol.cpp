/* Copyright (C) 2019-2023 Andrew D. Smith
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

[[maybe_unused]] static constexpr auto about = R"(
guessprotocol: a tool for guessing whether the protocol for generating
whole-genome bisulfite sequencing was wgbs, pbat or random pbat.
)";

#include "numerical_utils.hpp"

#include "OptionParser.hpp"

#include <bamxx.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iostream>
#include <iterator>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

static constexpr std::array<int, 96> nuc_to_idx = {
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 16
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 32
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 48
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 64
  4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,  // 80
  4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 96
};

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
    std::transform(std::cbegin(pr), std::cend(pr), std::begin(lpr),
                   [](const double x) { return std::log(x); });
  }

  auto
  operator()(const std::string &s) const -> double {
    static const auto acc = [&](const double x, const char c) {
      // NOLINTNEXTLINE(*-constant-array-index)
      const auto i = nuc_to_idx[static_cast<std::uint8_t>(c)];
      return i == 4 ? x : x + lpr[i];
    };
    return std::accumulate(std::cbegin(s), std::cend(s), 0.0, acc);
  };

  [[nodiscard]] auto
  tostring() const -> std::string {
    std::ostringstream oss;
    oss << "pr:\n";
    for (const auto i : pr)
      oss << i << '\n';
    oss << "log pr:\n";
    for (const auto i : lpr)
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
  std::uint64_t n_reads{};

  // wgbs_fraction is the probability that a read (for single-ended reads) or
  // the read1 of a read pair (for paired reads) is T-rich.
  double wgbs_fraction{};

  void
  evaluate() {
    const auto frac =
      static_cast<double>(n_reads_wgbs) / static_cast<double>(n_reads);

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

  [[nodiscard]] auto
  tostring() const -> std::string {
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
struct fastq_record {
  std::string name;
  std::string seq;
};

// see if two reads from two ends match to each other (they should
// have the same name)
static auto
mates(const std::size_t to_ignore_at_end,  // in case names have #0/1 name ends
      const fastq_record &a, const fastq_record &b) -> bool {
  assert(to_ignore_at_end < std::size(a.name));
  const auto name_end =
    std::cend(a.name) - to_ignore_at_end;  // NOLINT(*-narrowing-conversions)
  return std::equal(std::cbegin(a.name), name_end, std::cbegin(b.name));
}

// Read 4 lines one time from fastq and fill in the fastq_record structure
static auto
operator>>(bamxx::bgzf_file &s, fastq_record &r) -> bamxx::bgzf_file & {
  if (!getline(s, r.name))
    return s;

  if (r.name.empty() || r.name[0] != '@')
    throw std::runtime_error("failed to parse fastq name line");

  const auto nm_end = r.name.find_first_of(" \t");
  const auto nm_sz =
    (nm_end == std::string::npos ? std::size(r.name) : nm_end) - 1;
  r.name.erase(std::copy_n(std::cbegin(r.name) + 1, nm_sz, std::begin(r.name)),
               std::cend(r.name));

  if (!getline(s, r.seq))
    throw std::runtime_error("failed to parse fastq sequence line");

  std::string tmp;
  if (!getline(s, tmp))
    throw std::runtime_error("failed to parse fastq plus line");

  if (!getline(s, tmp))
    throw std::runtime_error("failed to parse fastq qual line");

  return s;
}

auto
main_guessprotocol(int argc, char *argv[]) -> int {  // NOLINT(*-avoid-c-arrays)
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
    static constexpr auto description =
      "guess bisulfite protocol for a library";

    bool verbose{};
    bool use_human{};
    std::string outfile;
    std::size_t reads_to_check{1000000};  // NOLINT(*-avoid-magic-numbers)
    std::size_t name_suffix_len{};
    double bisulfite_conversion_rate{0.98};  // NOLINT(*-avoid-magic-numbers)

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(argv[0],  // NOLINT(*-pointer-arithmetic)
                           description, "<end1-fastq> [<end2-fastq>]");
    opt_parse.add_opt("nreads", 'n', "number of reads in initial check", false,
                      reads_to_check);
    opt_parse.add_opt("ignore", 'i',
                      "length of read name suffix to ignore when matching",
                      false, name_suffix_len);
    opt_parse.add_opt("bisulfite", 'b', "bisulfite conversion rate", false,
                      bisulfite_conversion_rate);
    opt_parse.add_opt("human", 'H', "assume human genome", false, use_human);
    opt_parse.add_opt("output", 'o', "output file name", true, outfile);
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
    if (opt_parse.about_requested() || std::size(leftover_args) > 2) {
      std::cerr << opt_parse.about_message() << '\n';
      return EXIT_SUCCESS;
    }
    const std::vector<std::string> reads_files(leftover_args);
    /****************** END COMMAND LINE OPTIONS *****************/

    const auto base_comp = use_human ? human_base_comp : flat_base_comp;
    // if (use_human)
    //   base_comp = human_base_comp;

    nucleotide_model t_rich_model(base_comp, bisulfite_conversion_rate, true);
    nucleotide_model a_rich_model(base_comp, bisulfite_conversion_rate, false);

    guessprotocol_summary summary;
    summary.layout = std::size(reads_files) == 2 ? "paired" : "single";

    if (verbose) {
      if (std::size(reads_files) == 2)
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

    if (std::size(reads_files) == 2) {
      // input: paired-end reads with end1 and end2
      bamxx::bgzf_file in1(reads_files.front(), "r");
      if (!in1)
        throw std::runtime_error("cannot open file: " + reads_files.front());

      bamxx::bgzf_file in2(reads_files.back(), "r");
      if (!in2)
        throw std::runtime_error("cannot open file: " + reads_files.back());

      fastq_record r1, r2;
      while (in1 >> r1 && in2 >> r2 && summary.n_reads < reads_to_check) {
        summary.n_reads++;

        if (!mates(name_suffix_len, r1, r2))
          throw std::runtime_error("expected mates: " + r1.name + ", " +
                                   r2.name);

        const double ta = t_rich_model(r1.seq) + a_rich_model(r2.seq);
        const double at = a_rich_model(r1.seq) + t_rich_model(r2.seq);

        const auto prob_read1_t_rich = std::exp(ta - log_sum_log(ta, at));
        summary.n_reads_wgbs += prob_read1_t_rich;
      }
    }
    else {
      // input: single-end reads
      bamxx::bgzf_file in(reads_files.front(), "r");
      if (!in)
        throw std::runtime_error("cannot open file: " + reads_files.front());

      fastq_record r;
      while (in >> r && summary.n_reads < reads_to_check) {
        summary.n_reads++;

        const double t = t_rich_model(r.seq);
        const double a = a_rich_model(r.seq);

        const auto prob_t_rich = std::exp(t - log_sum_log(t, a));
        summary.n_reads_wgbs += prob_t_rich;
      }
    }

    summary.evaluate();

    std::ofstream out(outfile);
    if (!out)
      throw std::runtime_error("failed to open: " + outfile);
    out << summary.tostring() << '\n';
  }
  catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
