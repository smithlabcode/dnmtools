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

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>
#include <stdexcept>
#include <cmath>
#include <algorithm>

// smithlab headers
#include "OptionParser.hpp"
#include "smithlab_os.hpp"
#include "smithlab_utils.hpp"
#include "GenomicRegion.hpp"

using std::string;
using std::vector;
using std::cerr;
using std::cout;
using std::endl;
using std::istringstream;
using std::istream;
using std::ostream;
using std::ifstream;
using std::ofstream;
using std::runtime_error;

// Attemps to find the next significant CpG site. Returns true if one was found
// and flase otherwise.
static bool
read_next_significant_cpg(istream &cpg_stream, GenomicRegion &cpg,
                          double cutoff, bool &skipped_any, bool &n_sig_sites,
                          size_t &test_cov, size_t &test_meth,
                          size_t &rest_cov, size_t &rest_meth) {
  GenomicRegion region;
  skipped_any = false;
  n_sig_sites = false;
  string cpg_encoding;

  while (getline(cpg_stream, cpg_encoding)) {
    string record, chrom, name, sign;
    size_t position;
    double raw_pval, adjusted_pval, corrected_pval;

    istringstream iss(cpg_encoding);
    iss.exceptions(std::ios::failbit);
    iss >> chrom >> position >> sign >> name >> raw_pval
        >> adjusted_pval >> corrected_pval
        >> test_cov >> test_meth >> rest_cov >> rest_meth;

    if (corrected_pval >= 0.0 && corrected_pval < cutoff) {
      cpg.set_chrom(chrom);
      cpg.set_start(position);
      cpg.set_end(position + 1);
      n_sig_sites = (0 <= raw_pval && raw_pval < cutoff);
      return true;
    }
    skipped_any = true;
  }

  return false;
}

static void
merge(istream &cpg_stream, ostream &dmr_stream, double cutoff) {

  GenomicRegion dmr;
  dmr.set_name("dmr");

  size_t dmr_test_cov = 0;
  size_t dmr_test_meth = 0;
  size_t dmr_rest_cov = 0;
  size_t dmr_rest_meth = 0;

  size_t test_cov = 0;
  size_t test_meth = 0;
  size_t rest_cov = 0;
  size_t rest_meth = 0;

  // Find the first significant CpG, or terminate the function if none exist.
  bool skipped_last_cpg, n_sig_sites;
  if (!read_next_significant_cpg(cpg_stream, dmr, cutoff, skipped_last_cpg,
                                 n_sig_sites, test_cov, test_meth,
                                 rest_cov, rest_meth))
    return;

  dmr.set_score(n_sig_sites);
  dmr_test_cov += test_cov;
  dmr_test_meth += test_meth;
  dmr_rest_cov += rest_cov;
  dmr_rest_meth += rest_meth;

  GenomicRegion cpg;
  cpg.set_name("dmr");

  while (read_next_significant_cpg(cpg_stream, cpg, cutoff, skipped_last_cpg,
                                   n_sig_sites, test_cov, test_meth,
                                   rest_cov, rest_meth)) {

    if (skipped_last_cpg || cpg.get_chrom() != dmr.get_chrom()) {
      if (dmr.get_score() != 0)
        dmr_stream << dmr.get_chrom() << '\t'
                   << dmr.get_start() << '\t'
                   << dmr.get_end()   << '\t'
                   << dmr.get_name()  << '\t'
                   << dmr.get_score() << '\t'
                   << double(dmr_test_meth)/dmr_test_cov -
          double(dmr_rest_meth)/dmr_rest_cov << endl;
      dmr = cpg;
      dmr.set_score(n_sig_sites);
      dmr_test_cov = test_cov;
      dmr_test_meth = test_meth;
      dmr_rest_cov = rest_cov;
      dmr_rest_meth = rest_meth;
    }
    else {
      dmr.set_end(cpg.get_end());
      dmr.set_score(dmr.get_score() + n_sig_sites);
      dmr_test_cov += test_cov;
      dmr_test_meth += test_meth;
      dmr_rest_cov += rest_cov;
      dmr_rest_meth += rest_meth;
    }
  }
  if (dmr.get_score() != 0) {
    const double diff =
      static_cast<double>(dmr_test_meth)/dmr_test_cov -
      static_cast<double>(dmr_rest_meth)/dmr_rest_cov;
    dmr_stream << dmr.get_chrom() << '\t'
               << dmr.get_start() << '\t'
               << dmr.get_end()   << '\t'
               << dmr.get_name()  << '\t'
               << dmr.get_score() << '\t'
               << diff << endl;
  }
}

int
main_radmeth_merge(int argc, char *argv[]) {

  try {

    /* FILES */
    string outfile;
    double cutoff = 0.01;

    /**************** GET COMMAND LINE ARGUMENTS *************************/
    OptionParser opt_parse(strip_path(argv[0]),
                           "merge significantly differentially"
                           " methylated CpGs into DMRs",
                           "<bed-file-in-radmeth-format>");
    opt_parse.set_show_defaults();
    opt_parse.add_opt("output", 'o',
                      "output file (default: stdout)", false, outfile);
    opt_parse.add_opt("cutoff", 'p', "p-value cutoff", false , cutoff);
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl;
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
    const string bed_filename = leftover_args.front();
    /************************************************************************/

    ofstream of;
    if (!outfile.empty()) of.open(outfile);
    ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());

    ifstream in(bed_filename);
    if (!in)
      throw runtime_error("could not open file: " + bed_filename);

    merge(in, out, cutoff);

  }
  catch (const std::exception &e) {
    cerr << "ERROR: " << e.what() << endl;
    exit(EXIT_FAILURE);
  }
  return EXIT_SUCCESS;
}
