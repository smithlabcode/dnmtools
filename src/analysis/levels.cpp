/* levels: a program to compute coverage statistics, mutation rates,
 * and three different formulas for methylation levels described in
 * the paper:
 *
 *    'Leveling' the playing field for analyses of single-base
 *     resolution DNA methylomes
 *     Schultz, Schmitz & Ecker (TIG 2012)
 *
 * Note: the fractional methylation level calculated in this program
 * is inspired but different from the paper. What we are doing here is
 * using binomial test to determine significantly hyper/hypomethylated
 * sites, and only use these subset of sites to calculate methylation
 * level.
 *
 * Copyright (C) 2014-2023 University of Southern California and
 *                         Andrew D. Smith and Benjamin E Decato
 *
 * Authors: Andrew D. Smith and Benjamin E Decato
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

#include "LevelsCounter.hpp"
#include "MSite.hpp"
#include "OptionParser.hpp"
#include "bsutils.hpp"
#include "counts_header.hpp"
#include "smithlab_os.hpp"
#include "smithlab_utils.hpp"

#include <bamxx.hpp>

#include <cmath>
#include <exception>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

using std::cerr;
using std::endl;
using std::ios_base;
using std::runtime_error;
using std::string;
using std::vector;

using bamxx::bgzf_file;

enum class counts_file_format { ordinary, asym_cpg, sym_cpg };

static counts_file_format
guess_counts_file_format(const string &filename) {
  static const uint64_t n_lines_to_check = 10000;

  const bool has_counts_header = get_has_counts_header(filename);
  bgzf_file in(filename, "r");
  if (!in) throw ios_base::failure{"bad input file: " + filename};
  if (has_counts_header)
    skip_counts_header(in);

  bool found_non_cpg = false, found_asym_cpg = false;
  MSite curr_site, prev_site;
  uint64_t line_count = 0;

  while (read_site(in, curr_site) && (!found_non_cpg || !found_asym_cpg) &&
         line_count++ < n_lines_to_check) {
    found_non_cpg = found_non_cpg || !curr_site.is_cpg();
    found_asym_cpg =
      found_asym_cpg || (curr_site.is_cpg() && prev_site.is_cpg());
  }
  return (found_non_cpg ? counts_file_format::ordinary
                        : (found_asym_cpg ? counts_file_format::asym_cpg
                                          : counts_file_format::sym_cpg));
}

int
main_levels(int argc, const char **argv) {
  try {
    bool verbose = false;
    bool relaxed_mode = false;
    string outfile;

    const string description = "compute global summary of methylation levels";

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), description, "<counts-file>");
    opt_parse.add_opt("output", 'o', "output file (default: stdout)", false,
                      outfile);
    opt_parse.add_opt("alpha", 'a',
                      "alpha for confidence interval (default: 0.05)", false,
                      LevelsCounter::alpha);
    opt_parse.add_opt("relaxed", '\0',
                      "run on data that appears to have sites filtered", false,
                      relaxed_mode);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, verbose);
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (opt_parse.help_requested()) {
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
    const string meth_file = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/

    if (!is_msite_file(meth_file))
      throw runtime_error{"malformed counts file: " + meth_file};

    const counts_file_format guessed_format =
      guess_counts_file_format(meth_file);
    if (guessed_format != counts_file_format::ordinary) {
      if (verbose)
        cerr << "input might be only CpG sites ("
             << (guessed_format == counts_file_format::asym_cpg ? "not " : "")
             << "symmetric)" << endl;
      if (!relaxed_mode)
        throw runtime_error{
          "unexpected input format (consider using -relaxed)"};
    }

    const bool has_counts_header = get_has_counts_header(meth_file);
    bgzf_file in(meth_file, "r");
    if (!in) throw std::runtime_error("bad input file: " + meth_file);
    if (has_counts_header)
      skip_counts_header(in);

    LevelsCounter cpg("cpg");
    LevelsCounter cpg_symmetric("cpg_symmetric");
    LevelsCounter chh("chh");
    LevelsCounter cxg("cxg");
    LevelsCounter ccg("ccg");
    LevelsCounter cytosines("cytosines");

    MSite site, prev_site;

    while (read_site(in, site)) {
      if (site.chrom != prev_site.chrom)
        if (verbose) cerr << "processing " << site.chrom << endl;

      cytosines.update(site);

      if (site.is_cpg()) {
        cpg.update(site);
        if (site.is_mate_of(prev_site)) {
          site.add(prev_site);
          cpg_symmetric.update(site);
        }
      }
      else if (site.is_chh())
        chh.update(site);
      else if (site.is_ccg())
        ccg.update(site);
      else if (site.is_cxg())
        cxg.update(site);
      else
        throw runtime_error{"bad site context: " + site.context};

      prev_site = site;
    }

    std::ofstream of;
    if (!outfile.empty()) {
      of.open(outfile);
      if (!of) throw ios_base::failure{"bad output file: " + outfile};
    }
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

    out << cytosines << '\n'
        << cpg << '\n'
        << cpg_symmetric << '\n'
        << chh << '\n'
        << ccg << '\n'
        << cxg << '\n';
  }
  catch (const std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
