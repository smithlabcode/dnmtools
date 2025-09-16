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

#include "GenomicRegion.hpp"
#include "LevelsCounter.hpp"
#include "MSite.hpp"
#include "OptionParser.hpp"
#include "bsutils.hpp"
#include "smithlab_utils.hpp"
#include "xcounts_utils.hpp"

#include <bamxx.hpp>

#include <algorithm>
#include <charconv>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

using bamxx::bgzf_file;

namespace fs = std::filesystem;

static std::string
format_levels_counter(const LevelsCounter &lc) {
  // ...
  // (7) weighted mean methylation
  // (8) unweighted mean methylation
  // (9) fractional methylation
  // (10) number of sites in the region
  // (11) number of sites covered at least once
  // (12) number of observations in reads indicating methylation
  // (13) total number of observations from reads in the region
  std::ostringstream oss;
  // clang-format off
  oss << lc.mean_meth_weighted() << '\t'
      << lc.mean_meth() << '\t'
      << lc.fractional_meth() << '\t'
      << lc.total_sites << '\t'
      << lc.sites_covered << '\t'
      << lc.total_c << '\t'
      << (lc.total_c + lc.total_t) << '\t'
      << lc.total_meth << '\t'
      << lc.sites_covered;
  // clang-format on
  return oss.str();
}

static std::unordered_map<std::string, std::uint64_t>
get_chrom_sizes(const std::string &chrom_sizes_file) {
  std::unordered_map<std::string, std::uint64_t> chrom_sizes;

  std::ifstream in(chrom_sizes_file);
  if (!in)
    throw std::runtime_error("failed to open file: " + chrom_sizes_file);

  std::string line;
  while (getline(in, line)) {
    std::istringstream iss(line);
    std::string chrom_name{};
    std::uint64_t chrom_size{};
    if (!(iss >> chrom_name >> chrom_size))
      throw std::runtime_error("bad line in " + chrom_sizes_file + ":\n" +
                               line);
    std::string dummy;
    if (iss >> dummy)
      throw std::runtime_error("too many columns: " + line);

    if (chrom_sizes.find(chrom_name) != std::cend(chrom_sizes))
      throw std::runtime_error("repeated entry " + chrom_name + " in " +
                               chrom_sizes_file);

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
      throw std::runtime_error("bad line in " + chrom_sizes_file + ":\n" +
                               line);

    chrom_names.push_back(chrom_name);
  }
  return chrom_names;
}

static void
update(LevelsCounter &lc, const xcounts_entry &xce) {
  const std::uint64_t n_reads = xce.n_reads();
  if (n_reads > 0) {
    ++lc.sites_covered;
    lc.max_depth = std::max(lc.max_depth, n_reads);
    lc.total_c += xce.n_meth;
    lc.total_t += xce.n_unmeth;
    const auto meth = xce.frac();
    lc.total_meth += meth;
    double lower = 0.0, upper = 0.0;
    wilson_ci_for_binomial(lc.alpha, n_reads, meth, lower, upper);
    lc.called_meth += (lower > 0.5);
    lc.called_unmeth += (upper < 0.5);
  }
  ++lc.total_sites;
}

static void
process_chrom(const bool report_more_info, const char level_code,
              const std::string &chrom_name, const std::uint64_t chrom_size,
              const std::uint64_t bin_size,
              const std::vector<xcounts_entry> &sites, std::ostream &out) {
  GenomicRegion r(chrom_name, 0, 0, "CpG", 0.0, '+');

  std::uint64_t j = 0;
  for (auto i = 0ul; i < chrom_size; i += bin_size) {
    while (j < std::size(sites) && sites[j].pos < i)
      ++j;

    LevelsCounter lc;
    while (j < std::size(sites) && sites[j].pos < i + bin_size)
      update(lc, sites[j++]);

    r.set_start(i);
    r.set_end(std::min(i + bin_size, chrom_size));
    r.set_score(level_code == 'w' ? lc.mean_meth_weighted()
                                  : (level_code == 'u' ? lc.mean_meth()
                                                       : lc.fractional_meth()));
    r.set_name("CpG_" +
               std::to_string((level_code == 'w'
                                 ? lc.coverage()
                                 : (level_code == 'u' ? lc.sites_covered
                                                      : lc.total_called()))));
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
  GenomicRegion r(chrom_name, 0, 0, "CpG_0", 0.0, '+');
  LevelsCounter lc;
  const std::string lc_formatted = format_levels_counter(lc);
  for (auto i = 0ul; i < chrom_size; i += bin_size) {
    r.set_start(i);
    r.set_end(std::min(i + bin_size, chrom_size));
    out << r;
    if (report_more_info)
      out << '\t' << lc_formatted;
    out << '\n';
  }
}

int
main_cpgbins(int argc, char *argv[]) {
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

    static const std::string default_name_prefix = "X";

    bool verbose = false;
    bool report_more_info = false;
    std::uint32_t n_threads = 1;
    // std::uint64_t bin_size = 1000;
    size_t bin_size = 1000;  // ADS: for macOS gcc-14.2.0
    std::string level_code = "w";
    std::string outfile;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), description,
                           "<chrom-sizes> <xsym-file>");
    opt_parse.set_show_defaults();
    opt_parse.add_opt("output", 'o', "name of output file (default: stdout)",
                      false, outfile);
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

    if (!fs::is_regular_file(chrom_sizes_file))
      throw std::runtime_error("chromosome sizes file not a regular file: " +
                               chrom_sizes_file);

    if (!fs::is_regular_file(xcounts_file))
      throw std::runtime_error("xsym file not a regular file: " + xcounts_file);

    const auto sites_by_chrom = read_xcounts_by_chrom(n_threads, xcounts_file);
    const auto chrom_names = get_chrom_names(chrom_sizes_file);
    const auto chrom_sizes = get_chrom_sizes(chrom_sizes_file);

    std::ofstream of;
    if (!outfile.empty()) {
      of.open(outfile);
      if (!of)
        throw std::runtime_error("failed to open outfile: " + outfile);
    }
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

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
