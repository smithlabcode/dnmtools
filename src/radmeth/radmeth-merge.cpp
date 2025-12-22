/* Copyright (C) 2013-2023 University of Southern California and
 *                         Egor Dolzhenko
 *                         Andrew D Smith
 *
 * Authors: Andrew D. Smith and Egor Dolzhenko and Guilherme Sena
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
#include "OptionParser.hpp"

#include <cstdlib>
#include <exception>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// NOLINTBEGIN(*-narrowing-conversions)

// Attemps to find the next significant CpG site. Returns true if one was found
// and flase otherwise.
static bool
read_next_significant_cpg(std::istream &cpg_stream, Interval6 &cpg,
                          double cutoff, bool &skipped_any, bool &has_sig_sites,
                          std::size_t &test_cov, std::size_t &test_meth,
                          std::size_t &rest_cov, std::size_t &rest_meth) {
  Interval6 region;
  skipped_any = false;
  has_sig_sites = false;
  std::string line;

  while (std::getline(cpg_stream, line)) {
    std::string chrom, name, sign;
    std::size_t position{};
    double raw_pval{};
    double adjusted_pval{};
    double corrected_pval{};
    std::istringstream iss(line);
    if (!(iss >> chrom >> position >> sign >> name >> raw_pval >>
          adjusted_pval >> corrected_pval >> test_cov >> test_meth >>
          rest_cov >> rest_meth))
      throw std::runtime_error("failed to parse line:\n" + line);
    if (corrected_pval >= 0.0 && corrected_pval < cutoff) {
      cpg.chrom = chrom;
      cpg.start = position;
      cpg.stop = position + 1;
      has_sig_sites = (0.0 <= raw_pval && raw_pval < cutoff);
      return true;
    }
    skipped_any = true;
  }

  return false;
}

static void
merge(std::istream &cpg_stream, std::ostream &dmr_stream, double cutoff) {
  Interval6 dmr;
  dmr.name = "dmr";

  std::size_t dmr_test_cov{};
  std::size_t dmr_test_meth{};
  std::size_t dmr_rest_cov{};
  std::size_t dmr_rest_meth{};

  std::size_t test_cov{};
  std::size_t test_meth{};
  std::size_t rest_cov{};
  std::size_t rest_meth{};

  // Find the first significant CpG, or terminate the function if none exist.
  bool skipped_last_cpg{};
  bool has_sig_sites{};
  if (!read_next_significant_cpg(cpg_stream, dmr, cutoff, skipped_last_cpg,
                                 has_sig_sites, test_cov, test_meth, rest_cov,
                                 rest_meth))
    return;

  dmr.score = has_sig_sites ? 1.0 : 0.0;
  dmr_test_cov += test_cov;
  dmr_test_meth += test_meth;
  dmr_rest_cov += rest_cov;
  dmr_rest_meth += rest_meth;

  Interval6 cpg;
  cpg.name = "dmr";

  while (read_next_significant_cpg(cpg_stream, cpg, cutoff, skipped_last_cpg,
                                   has_sig_sites, test_cov, test_meth, rest_cov,
                                   rest_meth)) {
    if (skipped_last_cpg || cpg.chrom != dmr.chrom) {
      if (dmr.score != 0)
        dmr_stream << dmr.chrom << '\t' << dmr.start << '\t' << dmr.stop << '\t'
                   << dmr.name << '\t' << dmr.score << '\t'
                   << static_cast<double>(dmr_test_meth) / dmr_test_cov -
                        static_cast<double>(dmr_rest_meth) / dmr_rest_cov
                   << '\n';
      dmr = cpg;
      dmr.score += has_sig_sites ? 1.0 : 0.0;
      dmr_test_cov = test_cov;
      dmr_test_meth = test_meth;
      dmr_rest_cov = rest_cov;
      dmr_rest_meth = rest_meth;
    }
    else {
      dmr.stop = cpg.stop;
      dmr.score += has_sig_sites ? 1.0 : 0.0;
      dmr_test_cov += test_cov;
      dmr_test_meth += test_meth;
      dmr_rest_cov += rest_cov;
      dmr_rest_meth += rest_meth;
    }
  }
  if (dmr.score != 0) {
    const double diff = static_cast<double>(dmr_test_meth) / dmr_test_cov -
                        static_cast<double>(dmr_rest_meth) / dmr_rest_cov;
    dmr_stream << dmr.chrom << '\t' << dmr.start << '\t' << dmr.stop << '\t'
               << dmr.name << '\t' << dmr.score << '\t' << diff << '\n';
  }
}

int
main_radmeth_merge(int argc, char *argv[]) {  // NOLINT(*-avoid-c-arrays)
  try {
    std::string outfile;
    double cutoff = 0.01;  // NOLINT(*-avoid-magic-numbers)

    /**************** GET COMMAND LINE ARGUMENTS *************************/
    OptionParser opt_parse(argv[0],  // NOLINT(*-pointer-arithmetic)
                           "merge significantly differentially"
                           " methylated CpGs into DMRs",
                           "<bed-file-in-radmeth-format>");
    opt_parse.set_show_defaults();
    opt_parse.add_opt("output", 'o', "output file", true, outfile);
    opt_parse.add_opt("cutoff", 'p', "p-value cutoff", false, cutoff);
    std::vector<std::string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      std::cerr << opt_parse.help_message() << '\n';
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
    if (leftover_args.size() != 1) {
      std::cerr << opt_parse.help_message() << '\n';
      return EXIT_SUCCESS;
    }
    const std::string bed_filename = leftover_args.front();
    /************************************************************************/

    std::ofstream out(outfile);
    if (!out)
      throw std::runtime_error("failed to open: " + outfile);
    std::ifstream in(bed_filename);
    if (!in)
      throw std::runtime_error("failed to open: " + bed_filename);

    merge(in, out, cutoff);
  }
  catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

// NOLINTEND(*-narrowing-conversions)
