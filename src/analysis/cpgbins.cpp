/* cpgbins: average CpG methylation in genomic bins
 *
 * Copyright (C) 2024 Andrew D. Smith
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

#include "Interval6.hpp"
#include "LevelsCounter.hpp"
#include "OptionParser.hpp"
#include "bsutils.hpp"
#include "xcounts_utils.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <stdlib.h>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

static std::unordered_map<std::string, std::uint64_t>
get_chrom_sizes(const std::string &chrom_sizes_file) {
  std::ifstream in(chrom_sizes_file);
  if (!in)
    throw std::runtime_error("failed to open file: " + chrom_sizes_file);

  std::unordered_map<std::string, std::uint64_t> chrom_sizes;

  std::string line;
  while (getline(in, line)) {
    std::istringstream iss(line);
    std::string chrom_name{};
    std::uint64_t chrom_size{};
    if (!(iss >> chrom_name >> chrom_size))
      throw std::runtime_error("bad chromosome sizes line:\n" + line);
    std::string dummy;
    if (iss >> dummy)
      throw std::runtime_error("too many columns: " + line);

    if (chrom_sizes.find(chrom_name) != std::cend(chrom_sizes))
      throw std::runtime_error("repeated entry in chromosome sizes: " + line);

    chrom_sizes[chrom_name] = chrom_size;
  }
  return chrom_sizes;
}

static std::vector<std::string>
get_chrom_names(const std::string &chrom_sizes_file) {
  std::ifstream in(chrom_sizes_file);
  if (!in)
    throw std::runtime_error("failed to open file: " + chrom_sizes_file);

  std::vector<std::string> chrom_names;

  std::string line;
  while (getline(in, line)) {
    std::istringstream iss(line);
    std::string chrom_name{};
    if (!(iss >> chrom_name))
      throw std::runtime_error("bad line in chromosome sizes file: " + line);
    chrom_names.push_back(chrom_name);
  }
  return chrom_names;
}

static void
update(LevelsCounter &lc, const xcounts_entry &xce) {
  static constexpr auto one_half = 0.5;
  const std::uint64_t n_reads = xce.n_reads();
  if (n_reads > 0) {
    ++lc.sites_covered;
    lc.max_depth = std::max(lc.max_depth, n_reads);
    lc.total_c += xce.n_meth;
    lc.total_t += xce.n_unmeth;
    const auto meth = xce.frac();
    lc.total_meth += meth;
    double lower{}, upper{};
    wilson_ci_for_binomial(lc.alpha, static_cast<double>(n_reads), meth, lower,
                           upper);
    lc.called_meth += (lower > one_half);
    lc.called_unmeth += (upper < one_half);
  }
  ++lc.total_sites;
}

[[nodiscard]] static std::string
get_name(const bool report_more_info, const char level_code,
         const LevelsCounter &lc) {
  static const std::string tag = "CpG";
  return tag + (report_more_info
                  ? ""
                  : "_" + std::to_string(level_code == 'w'
                                           ? lc.coverage()
                                           : (level_code == 'u'
                                                ? lc.sites_covered
                                                : lc.total_called())));
}

static void
process_chrom(const bool report_more_info, const char level_code,
              const std::string &chrom_name, const std::uint64_t chrom_size,
              const std::uint64_t bin_size,
              const std::vector<xcounts_entry> &sites, std::ostream &out) {
  Interval6 r(chrom_name, 0, 0, "CpG", 0.0, '+');
  std::uint64_t j{};
  for (auto i = 0ul; i < chrom_size; i += bin_size) {
    while (j < std::size(sites) && sites[j].pos < i)
      ++j;

    LevelsCounter lc;
    while (j < std::size(sites) && sites[j].pos < i + bin_size)
      update(lc, sites[j++]);

    r.start = i;
    r.stop = std::min(i + bin_size, chrom_size);
    r.score = level_code == 'w'
                ? lc.mean_meth_weighted()
                : (level_code == 'u' ? lc.mean_meth() : lc.fractional_meth());

    r.name = get_name(report_more_info, level_code, lc);

    out << r;
    if (report_more_info)
      out << '\t' << format_levels_counter(lc);
    out << '\n';
  }
}

static void
process_chrom(const bool report_more_info, const std::string &chrom_name,
              const std::uint64_t chrom_size, const std::uint64_t bin_size,
              std::ostream &out) {
  Interval6 r(chrom_name, 0, 0, "CpG_0", 0.0, '+');
  LevelsCounter lc;
  const std::string lc_formatted = format_levels_counter(lc);
  for (auto i = 0ul; i < chrom_size; i += bin_size) {
    r.start = i;
    r.stop = std::min(i + bin_size, chrom_size);
    out << r;
    if (report_more_info)
      out << '\t' << lc_formatted;
    out << '\n';
  }
}

int
main_cpgbins(int argc, char *argv[]) {  // NOLINT(*-avoid-c-arrays)
  try {
    static const std::string description = R"""(
Compute average site methylation levels in each non-overlapping
genomic bin of the specified size. The 5th column (the "score" column
in BED format) is determined by the '-l' or '-level' argument.

Columns (beyond the first 6) in the BED format output:
(7) weighted mean methylation
(8) unweighted mean methylation
(9) fractional methylation
(10) number of sites in the region
(11) number of sites covered at least once
(12) number of observations in reads indicating methylation
(13) total number of observations from reads in the region
)""";

    bool verbose = false;
    bool report_more_info = false;
    std::int32_t n_threads = 1;
    // std::uint64_t bin_size = 1000;

    // ADS: for macOS gcc-14.2.0
    std::size_t bin_size = 1000;  // NOLINT(*-avoid-magic-numbers)
    std::string level_code = "w";
    std::string outfile;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(argv[0],  // NOLINT(*-pointer-arithmetic)
                           description, "<chrom-sizes> <xsym-file>");
    opt_parse.set_show_defaults();
    opt_parse.add_opt("output", 'o', "name of output file (default: stdout)",
                      true, outfile);
    opt_parse.add_opt("bin", 'b', "bin size in base pairs", false, bin_size);
    opt_parse.add_opt("threads", 't', "number of threads", false, n_threads);
    opt_parse.add_opt("level", 'l',
                      "the level to report as score column "
                      "in bed format output (w, u or f)",
                      false, level_code);
    opt_parse.add_opt("more-levels", 'M', "report more methylation information",
                      false, report_more_info);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, verbose);
    std::vector<std::string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested() ||
        opt_parse.about_requested()) {
      std::cerr << opt_parse.help_message() << opt_parse.about_message_raw()
                << '\n';
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      std::cerr << opt_parse.option_missing_message() << '\n';
      return EXIT_FAILURE;
    }
    if (std::size(leftover_args) != 2) {
      std::cerr << opt_parse.help_message() << '\n';
      return EXIT_FAILURE;
    }
    if (level_code != "w" && level_code != "u" && level_code != "f") {
      std::cerr << "selected level must be in {w, u, f}" << '\n';
      return EXIT_FAILURE;
    }
    const std::string chrom_sizes_file = leftover_args.front();
    const std::string xcounts_file = leftover_args.back();
    /****************** END COMMAND LINE OPTIONS *****************/

    if (n_threads <= 0) {
      std::cerr << "number of threads must be at least 1\n";
      return EXIT_FAILURE;
    }

    if (!std::filesystem::is_regular_file(chrom_sizes_file))
      throw std::runtime_error("chromosome sizes file not a regular file: " +
                               chrom_sizes_file);

    if (!std::filesystem::is_regular_file(xcounts_file))
      throw std::runtime_error("xsym file not a regular file: " + xcounts_file);

    const auto sites_by_chrom = read_xcounts_by_chrom(n_threads, xcounts_file);
    const auto chrom_names = get_chrom_names(chrom_sizes_file);
    const auto chrom_sizes = get_chrom_sizes(chrom_sizes_file);

    std::ofstream out(outfile);
    if (!out)
      throw std::runtime_error("failed to open outfile: " + outfile);

    for (const auto &chrom_name : chrom_names) {
      const auto sites = sites_by_chrom.find(chrom_name);
      const auto chrom_size = chrom_sizes.find(chrom_name);
      if (chrom_size == std::cend(chrom_sizes))
        throw std::runtime_error("failed chrom size lookup");
      if (sites != std::cend(sites_by_chrom))
        process_chrom(report_more_info, level_code[0], chrom_name,
                      chrom_size->second, bin_size, sites->second, out);
      else
        process_chrom(report_more_info, chrom_name, chrom_size->second,
                      bin_size, out);
    }
  }
  catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
