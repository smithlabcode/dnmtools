/* Copyright (C) 2013-2023 University of Southern California and
 *                         Egor Dolzhenko
 *                         Andrew D Smith
 *                         Guilherme Sena
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

#include <gsl/gsl_cdf.h> // GSL header for chisqrd distribution

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// smithlab headers
#include "GenomicRegion.hpp"
#include "OptionParser.hpp"
#include "smithlab_os.hpp"
#include "smithlab_utils.hpp"

#include "radmeth_model.hpp"
#include "radmeth_optimize.hpp"

using std::begin;
using std::cerr;
using std::cout;
using std::distance;
using std::end;
using std::endl;
using std::istringstream;
using std::runtime_error;
using std::sort;
using std::string;
using std::vector;

static std::istream &
operator>>(std::istream &is, Design &design) {
  string header_encoding;
  getline(is, header_encoding);

  istringstream header_is(header_encoding);
  string header_name;
  while (header_is >> header_name) design.factor_names.push_back(header_name);

  string row;
  while (getline(is, row)) {
    if (row.empty()) continue;

    istringstream row_is(row);
    string token;
    row_is >> token;
    design.sample_names.push_back(token);

    vector<double> matrix_row;
    while (row_is >> token) {
      if (token.length() == 1 && (token == "0" || token == "1"))
        matrix_row.push_back(token == "1");
      else
        throw runtime_error("only binary factor levels are allowed:\n" + row);
    }

    if (matrix_row.size() != design.n_factors())
      throw runtime_error(
        "each row must have as many columns as factors:\n" +
        row);

    design.matrix.push_back(vector<double>());
    swap(design.matrix.back(), matrix_row);
  }
  return is;
}

static std::ostream &
operator<<(std::ostream &out, const Design &design) {
  for (size_t factor = 0; factor < design.factor_names.size(); ++factor) {
    out << design.factor_names[factor];
    if (factor + 1 != design.factor_names.size()) out << "\t";
  }
  out << endl;

  for (size_t i = 0; i < design.n_samples(); ++i) {
    out << design.sample_names[i];
    for (size_t factor = 0; factor < design.factor_names.size(); ++factor)
      out << "\t" << design.matrix[i][factor];
    out << "\n";
  }
  return out;
}

static void
remove_factor(Design &design, const size_t factor) {
  design.factor_names.erase(begin(design.factor_names) + factor);
  for (size_t i = 0; i < design.n_samples(); ++i)
    design.matrix[i].erase(begin(design.matrix[i]) + factor);
}

// Parses a natural number from its string representation. Throws exception if
// the string does not encode one.
static size_t
parse_natural_number(string encoding) {
  istringstream iss(encoding);
  size_t number;
  iss >> number;
  if (!iss)
    throw runtime_error("The token \"" + encoding +
                        "\" does not encode a natural number");
  return number;
}

static std::istream &
operator>>(std::istream &table_encoding, SiteProportions &props) {
  props.chrom.clear();
  props.position = 0;
  props.strand.clear();
  props.context.clear();
  props.meth.clear();
  props.total.clear();

  string row_encoding;
  getline(table_encoding, row_encoding);

  // Skip lines contining only the newline character (e.g. the last line of the
  // proportion table).
  if (row_encoding.empty()) return table_encoding;

  // Get the row name (which must be specified like this: "chr:position") and
  // parse it.
  istringstream row_stream(row_encoding);
  string row_name_encoding;
  row_stream >> row_name_encoding;

  // Every row must start an identifier consisiting of genomic loci of the
  // corresponding site. Here we check this identifier has the correct number
  // of colons.
  const size_t num_colon =
    count(row_name_encoding.begin(), row_name_encoding.end(), ':');

  if (num_colon != 3)
    throw runtime_error(
      "Each row in the count table must start with "
      "a line chromosome:position:strand:context. Got \"" +
      row_name_encoding + "\" instead.");

  // First parse the row identifier.
  istringstream name_stream(row_name_encoding);
  getline(name_stream, props.chrom, ':');

  if (props.chrom.empty())
    throw runtime_error("Error parsing " + row_name_encoding +
                        ": chromosome name is missing.");

  string position_encoding;

  getline(name_stream, position_encoding, ':');
  props.position = parse_natural_number(position_encoding);
  getline(name_stream, props.strand, ':');
  getline(name_stream, props.context, ':');

  // After parsing the row identifier, parse count proportions.
  size_t total_count, meth_count;

  while (row_stream >> total_count >> meth_count) {
    props.total.push_back(total_count);
    props.meth.push_back(meth_count);
  }

  if (!row_stream.eof())
    throw runtime_error("Some row entries are not natural numbers: " +
                        row_stream.str());

  if (props.total.size() != props.meth.size())
    throw runtime_error(
      "This row does not encode proportions"
      "correctly:\n" +
      row_encoding);
  return table_encoding;
}

static bool
fit_regression_model(Regression &r) {
  vector<double> initial_params;
  return fit_regression_model(r, initial_params);
}

/***************** RADMETH ALGORITHM *****************/
static bool
consistent_sample_names(const Regression &reg, const string &header) {
  istringstream iss(header);
  auto nm_itr(begin(reg.design.sample_names));
  const auto nm_end(end(reg.design.sample_names));
  string token;
  while (iss >> token && nm_itr != nm_end)
    if (token != *nm_itr++) return false;
  return true;
}

// Given the maximum likelihood estimates of the full and reduced models, the
// function outputs the p-value of the log-likelihood ratio. *Note* that it is
// assumed that the reduced model has one fewer factor than the reduced model.
static double
loglikratio_test(double null_loglik, double full_loglik) {
  // The log-likelihood ratio statistic.
  const double log_lik_stat = -2 * (null_loglik - full_loglik);

  // It is assumed that null model has one fewer factor than the full
  // model. Hence the number of degrees of freedom is 1.
  const size_t degrees_of_freedom = 1;

  // Log-likelihood ratio statistic has a chi-sqare distribution.
  const double chisq_p = gsl_cdf_chisq_P(log_lik_stat, degrees_of_freedom);
  const double p_value = 1.0 - chisq_p;

  return p_value;
}

static bool
has_low_coverage(const Regression &reg, const size_t test_fac) {
  bool cvrd_in_test_fact_smpls = false;
  for (size_t i = 0; i < reg.n_samples() && !cvrd_in_test_fact_smpls; ++i)
    cvrd_in_test_fact_smpls =
      (reg.design.matrix[i][test_fac] == 1 && reg.props.total[i] != 0);

  bool cvrd_in_other_smpls = false;
  for (size_t i = 0; i < reg.n_samples() && !cvrd_in_other_smpls; ++i)
    cvrd_in_other_smpls =
      (reg.design.matrix[i][test_fac] != 1 && reg.props.total[i] != 0);

  return !cvrd_in_test_fact_smpls || !cvrd_in_other_smpls;
}

static bool
has_extreme_counts(const Regression &reg) {
  bool is_maximally_methylated = true;
  for (size_t i = 0; i < reg.n_samples() && is_maximally_methylated; ++i)
    is_maximally_methylated = (reg.props.total[i] == reg.props.meth[i]);

  bool is_unmethylated = true;
  for (size_t i = 0; i < reg.n_samples() && is_unmethylated; ++i)
    is_unmethylated = (reg.props.meth[i] == 0.0);

  return is_maximally_methylated || is_unmethylated;
}

/***********************************************************************
 * Run beta-binoimial regression using the specified table with
 * proportions and design matrix
 */
int
main_radmeth(int argc, const char **argv) {
  try {
    static const string description =
      "calculate differential methylation scores";

    string outfile;
    string test_factor_name;
    bool VERBOSE = false;
    bool more_na_info = false;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), description,
                           "<design-matrix> <data-matrix>");
    opt_parse.set_show_defaults();
    opt_parse.add_opt("out", 'o', "output file (default: stdout)", false,
                      outfile);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    opt_parse.add_opt(
      "na-info", 'n',
      "if a p-value is not calculated, print NAs in more detail: "
      "low count (NA_LOW_COV) extreme values (NA_EXTREME_CNT) or "
      "numerical errors in likelihood ratios (NA)",
      false, more_na_info);
    opt_parse.add_opt("factor", 'f', "a factor to test", true,
                      test_factor_name);

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
    if (leftover_args.size() != 2) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string design_filename(leftover_args.front());
    const string table_filename(leftover_args.back());
    /****************** END COMMAND LINE OPTIONS *****************/

    if (VERBOSE) cerr << "design table filename: " << design_filename << endl;
    std::ifstream design_file(design_filename);
    if (!design_file)
      throw runtime_error("could not open file: " + design_filename);

    // initialize full design matrix from file
    Regression full_regression;
    design_file >> full_regression.design;

    if (VERBOSE) cerr << full_regression.design << endl;

    // Check that provided test factor name exists and find its index.
    // Identify test factors with their indexes to simplify naming
    auto test_fact_it =
      find(begin(full_regression.design.factor_names),
           end(full_regression.design.factor_names), test_factor_name);

    if (test_fact_it == end(full_regression.design.factor_names))
      throw runtime_error(test_factor_name +
                          " factor not part of design specification");

    const size_t test_factor =
      distance(begin(full_regression.design.factor_names), test_fact_it);

    Regression null_regression;
    null_regression.design = full_regression.design;
    remove_factor(null_regression.design, test_factor);
    // ADS: done setup for the model


    // ADS: open the data table file
    std::ifstream table_file(table_filename);
    if (!table_file)
      throw runtime_error("could not open file: " + table_filename);

    // Make sure that the first line of the proportion table file contains
    // names of the samples. Throw an exception if the names or their order
    // in the proportion table does not match those in the full design matrix.
    string sample_names_header;
    getline(table_file, sample_names_header);

    if (!consistent_sample_names(full_regression, sample_names_header))
      // .design.sample_names != split_whitespace(sample_names_header))
      throw runtime_error("header:\n" + sample_names_header + "\n" +
                          "does not match factor names or their order in the\n"
                          "design matrix. Check that the design matrix and\n"
                          "the proportion table are correctly formatted.");

    const size_t n_samples = full_regression.design.n_samples();

    // ADS: now open the output file because each row is processed
    // sequentially
    std::ofstream of;
    if (!outfile.empty()) of.open(outfile);
    std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());

    // Performing the log-likelihood ratio test on proportions from each row
    // of the proportion table.
    while (table_file >> full_regression.props) {
      if (full_regression.props.total.size() != n_samples)
        throw runtime_error("found row with wrong number of columns");

      size_t coverage_factor = 0, coverage_rest = 0, meth_factor = 0,
             meth_rest = 0;

      for (size_t s = 0; s < n_samples; ++s) {
        if (full_regression.design.matrix[s][test_factor] != 0) {
          coverage_factor += full_regression.props.total[s];
          meth_factor += full_regression.props.meth[s];
        }
        else {
          coverage_rest += full_regression.props.total[s];
          meth_rest += full_regression.props.meth[s];
        }
      }

      out << full_regression.props.chrom << "\t"
          << full_regression.props.position << "\t"
          << full_regression.props.strand << "\t"
          << full_regression.props.context << "\t";

      // Do not perform the test if there's no coverage in either all
      // case or all control samples. Also do not test if the site is
      // completely methylated or completely unmethylated across all
      // samples.
      if (has_low_coverage(full_regression, test_factor)) {
        out << ((more_na_info) ? "NA_LOW_COV" : "NA");
      }
      else if (has_extreme_counts(full_regression)) {
        out << ((more_na_info) ? "NA_EXTREME_CNT" : "NA");
      }
      else {
        fit_regression_model(full_regression);
        null_regression.props = full_regression.props;
        fit_regression_model(null_regression);
        const double p_value = loglikratio_test(null_regression.max_loglik,
                                                full_regression.max_loglik);

        // If error occured in fitting (p-val = nan or -nan).
        if (p_value != p_value)
          out << "NA";
        else
          out << p_value;
      }
      out << "\t" << coverage_factor << "\t" << meth_factor << "\t"
          << coverage_rest << "\t" << meth_rest << endl;
    }
  }
  catch (const std::exception &e) {
    cerr << "ERROR: " << e.what() << endl;
    exit(EXIT_FAILURE);
  }
  return EXIT_SUCCESS;
}
