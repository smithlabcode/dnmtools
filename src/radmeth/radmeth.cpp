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

#include <gsl/gsl_cdf.h>  // GSL header for chisqrd distribution

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

// smithlab headers
#include "GenomicRegion.hpp"
#include "OptionParser.hpp"
#include "smithlab_os.hpp"
#include "smithlab_utils.hpp"

#include "radmeth_model.hpp"
#include "radmeth_optimize.hpp"

using std::begin;
using std::cout;
using std::distance;
using std::end;

struct file_progress {
  double one_hundred_over_filesize{};
  std::size_t prev_offset{};
  explicit file_progress(const std::string &filename) :
    one_hundred_over_filesize{100.0 / std::filesystem::file_size(filename)} {}
  void
  operator()(std::ifstream &in) {
    const std::size_t curr_offset = in.tellg() * one_hundred_over_filesize;
    if (curr_offset <= prev_offset)
      return;
    std::cerr << "[progress: " << std::setw(3) << curr_offset
              << (curr_offset == 100 ? "%]\n" : "%]\r");
    prev_offset = (curr_offset == 100) ? std::numeric_limits<std::size_t>::max()
                                       : curr_offset;
  }
};

static void
remove_factor(Design &design, const std::size_t factor_idx) {
  design.factor_names.erase(std::begin(design.factor_names) + factor_idx);
  for (std::size_t i = 0; i < design.n_samples(); ++i)
    design.matrix[i].erase(std::begin(design.matrix[i]) + factor_idx);
}

/***************** RADMETH ALGORITHM *****************/
static bool
consistent_sample_names(const Regression &reg, const std::string &header) {
  std::istringstream iss(header);
  auto nm_itr(std::begin(reg.design.sample_names));
  const auto nm_end(std::end(reg.design.sample_names));
  std::string token;
  while (iss >> token && nm_itr != nm_end)
    if (token != *nm_itr++)
      return false;
  return true;
}

// Given the maximum likelihood estimates of the full and reduced models, the
// function outputs the p-value of the log-likelihood ratio. *Note* that it is
// assumed that the reduced model has one fewer factor than the reduced model.
static double
llr_test(double null_loglik, double full_loglik) {
  // The log-likelihood ratio statistic.
  const double log_lik_stat = -2 * (null_loglik - full_loglik);

  // It is assumed that null model has one fewer factor than the full
  // model. Hence the number of degrees of freedom is 1.
  const std::size_t degrees_of_freedom = 1;

  // Log-likelihood ratio statistic has a chi-sqare distribution.
  const double chisq_p = gsl_cdf_chisq_P(log_lik_stat, degrees_of_freedom);
  const double p_value = 1.0 - chisq_p;

  return p_value;
}

static bool
has_low_coverage(const Regression &reg, const std::size_t test_fac) {
  bool cvrd_in_test_fact_smpls = false;
  for (std::size_t i = 0; i < reg.n_samples() && !cvrd_in_test_fact_smpls; ++i)
    cvrd_in_test_fact_smpls =
      (reg.design.matrix[i][test_fac] == 1 && reg.props.mc[i].n_reads != 0);

  bool cvrd_in_other_smpls = false;
  for (std::size_t i = 0; i < reg.n_samples() && !cvrd_in_other_smpls; ++i)
    cvrd_in_other_smpls =
      (reg.design.matrix[i][test_fac] != 1 && reg.props.mc[i].n_reads != 0);

  return !cvrd_in_test_fact_smpls || !cvrd_in_other_smpls;
}

[[nodiscard]] static bool
has_extreme_counts(const Regression &reg) {
  const auto &mc = reg.props.mc;

  bool full_meth = true;
  for (auto i = std::cbegin(mc); i != std::cend(mc) && full_meth; ++i)
    full_meth = (i->n_reads == i->n_meth);

  bool no_meth = true;
  for (auto i = std::cbegin(mc); i != std::cend(mc) && no_meth; ++i)
    no_meth = (i->n_meth == 0.0);

  return full_meth || no_meth;
}

[[nodiscard]] static bool
has_two_values(const Regression &reg, const std::size_t test_factor) {
  const auto &v = reg.design.tmatrix[test_factor];
  for (auto i = std::cbegin(v); i != std::cend(v); ++i)
    if (*i != v[0])
      return true;
  return false;
}

[[nodiscard]] static std::uint32_t
get_test_factor_idx(const Regression &model, const std::string &test_factor) {
  const auto &factors = model.design.factor_names;
  const auto itr =
    std::find(std::cbegin(factors), std::cend(factors), test_factor);

  if (itr == std::cend(factors))
    throw std::runtime_error("Factor not part of design: " + test_factor);

  return std::distance(std::cbegin(factors), itr);
}

static void
read_design(const bool verbose, const std::string &design_filename,
            Design &design) {
  if (verbose)
    std::cerr << "design table filename: " << design_filename << std::endl;
  std::ifstream design_file(design_filename);
  if (!design_file)
    throw std::runtime_error("could not open file: " + design_filename);

  // initialize full design matrix from file
  design_file >> design;

  if (verbose)
    std::cerr << design << std::endl;
}

enum class row_status : std::uint8_t {
  ok,
  na,
  na_low_cov,
  na_extreme_cnt,
};

[[nodiscard]] static std::vector<double>
drop_idx(const std::vector<double> &v, const std::size_t to_drop) {
  std::vector<double> u;
  u.reserve(std::size(v) - 1);
  for (auto i = 0u; i < std::size(v); ++i)
    if (i != to_drop)
      u.push_back(v[i]);
  return u;
}

static inline void
get_chunk_bounds(const std::uint32_t n_elements, const std::uint32_t n_chunks,
                 std::vector<std::pair<std::uint32_t, std::uint32_t>> &chunks) {
  const std::uint32_t q = n_elements / n_chunks;
  const std::uint32_t r = n_elements - q * n_chunks;
  std::uint32_t block_start{};
  for (std::size_t i = 0; i < n_chunks; ++i) {
    const auto sz = (i < r) ? q + 1 : q;
    const auto block_end = block_start + sz;
    chunks[i] = {block_start, block_end};
    block_start = block_end;
  }
}

static void
radmeth(const bool show_progress, const bool more_na_info,
        const std::uint32_t n_threads, const std::string &table_filename,
        const std::string &outfile, const Regression &alt_model,
        const Regression &null_model, const std::uint32_t test_factor_idx) {
  static constexpr auto prefix_fmt = "%s\t%ld\t%c\t%s\t";
  static constexpr auto suffix_fmt = "\t%ld\t%ld\t%ld\t%ld\n";
  static constexpr auto buf_size = 1024;
  static constexpr auto max_lines = 16384;

  // ADS: open the data table file
  std::ifstream table_file(table_filename);
  if (!table_file)
    throw std::runtime_error("could not open file: " + table_filename);

  // Make sure that the first line of the proportion table file contains
  // names of the samples. Throw an exception if the names or their order
  // in the proportion table does not match those in the full design matrix.
  std::string sample_names_header;
  std::getline(table_file, sample_names_header);

  if (!consistent_sample_names(alt_model, sample_names_header))
    throw std::runtime_error(
      "header:\n" + sample_names_header + "\n" +
      "does not match factor names or their order in the\n"
      "design matrix. Check that the design matrix and\n"
      "the proportion table are correctly formatted.");

  const std::size_t n_samples = alt_model.design.n_samples();

  std::ofstream out(outfile);
  if (!out)
    throw std::runtime_error("failed to open output file: " + outfile);

  file_progress progress{table_filename};

  std::vector<std::vector<char>> bufs(max_lines,
                                      std::vector<char>(buf_size, 0));
  std::vector<int> n_bytes(max_lines, 0);
  std::vector<std::string> lines(max_lines);

  std::vector<Regression> alt_models(n_threads, alt_model);
  std::vector<Regression> null_models(n_threads, null_model);

  std::vector<std::pair<std::uint32_t, std::uint32_t>> chunks(n_threads);

  // Iterate over rows in the file and do llr test on proportions from
  // each. Do this in sets of rows to avoid having to spawn too many threads.
  while (true) {
    std::uint32_t n_lines = 0;
    while (n_lines < max_lines && std::getline(table_file, lines[n_lines]))
      ++n_lines;
    if (n_lines == 0)
      break;

    get_chunk_bounds(n_lines, n_threads, chunks);

    std::vector<std::thread> threads;
    for (auto thread_id = 0u; thread_id < n_threads; ++thread_id) {
      threads.emplace_back([&, thread_id] {
        const auto &[chunk_beg, chunk_end] = chunks[thread_id];
        auto &t_alt_model = alt_models[thread_id];
        auto &t_null_model = null_models[thread_id];
        for (auto b = chunk_beg; b != chunk_end; ++b) {
          t_alt_model.props.parse(lines[b]);
          if (t_alt_model.props_size() != n_samples)
            throw std::runtime_error("found row with wrong number of columns");

          const auto [p_val, status] = [&]() -> std::tuple<double, row_status> {
            // Skip the test if (1) no coverage in all cases or in all controls,
            // or (2) the site is completely methylated or completely
            // unmethylated across all samples.
            if (has_low_coverage(t_alt_model, test_factor_idx))
              return std::tuple{1.0, row_status::na_low_cov};

            if (has_extreme_counts(t_alt_model))
              return std::tuple{1.0, row_status::na_extreme_cnt};

            std::vector<double> alternate_params;
            fit_regression_model(t_alt_model, alternate_params);
            t_null_model.props = t_alt_model.props;

            auto null_params = drop_idx(alternate_params, test_factor_idx);

            fit_regression_model(t_null_model, null_params);
            const double p_value =
              llr_test(t_null_model.max_loglik, t_alt_model.max_loglik);

            return (p_value != p_value) ? std::tuple{1.0, row_status::na}
                                        : std::tuple{p_value, row_status::ok};
          }();

          std::size_t n_reads_factor = 0;
          std::size_t n_reads_others = 0;
          std::size_t n_meth_factor = 0;
          std::size_t n_meth_others = 0;

          const auto &mc = t_alt_model.props.mc;
          const auto &vec = t_alt_model.design.tmatrix[test_factor_idx];
          for (std::size_t s = 0; s < n_samples; ++s)
            if (vec[s] != 0) {
              n_reads_factor += mc[s].n_reads;
              n_meth_factor += mc[s].n_meth;
            }
            else {
              n_reads_others += mc[s].n_reads;
              n_meth_others += mc[s].n_meth;
            }

          n_bytes[b] = [&] {
            // clang-format off
            const int n_prefix_bytes = std::sprintf(bufs[b].data(), prefix_fmt,
                                                    t_alt_model.props.chrom.data(),
                                                    t_alt_model.props.position,
                                                    t_alt_model.props.strand,
                                                    t_alt_model.props.context.data());
            // clang-format on
            if (n_prefix_bytes < 0)
              return n_prefix_bytes;

            auto cursor = bufs[b].data() + n_prefix_bytes;

            const int n_pval_bytes = [&] {
              if (status == row_status::ok)
                return std::sprintf(cursor, "%.6g", p_val);
              if (!more_na_info || status == row_status::na)
                return std::sprintf(cursor, "NA");
              if (status == row_status::na_extreme_cnt)
                return std::sprintf(cursor, "NA_EXTREME_CNT");
              // if (status == row_status::na_low_cov)
              return std::sprintf(cursor, "NA_LOW_COV");
            }();

            if (n_pval_bytes < 0)
              return n_pval_bytes;

            cursor += n_pval_bytes;

            // clang-format off
          const int n_suffix_bytes = std::sprintf(cursor, suffix_fmt,
                                                  n_reads_factor,
                                                  n_meth_factor,
                                                  n_reads_others,
                                                  n_meth_others);
            // clang-format on
            if (n_suffix_bytes < 0)
              return n_suffix_bytes;

            return n_prefix_bytes + n_pval_bytes + n_suffix_bytes;
          }();
        }
      });
    }

    for (auto &thread : threads)
      thread.join();

    for (auto i = 0u; i < n_lines; ++i) {
      if (n_bytes[i] < 0)
        throw std::runtime_error("failed to write to output buffer");
      out.write(bufs[i].data(), n_bytes[i]);
    }

    if (show_progress)
      progress(table_file);

    if (n_lines < n_threads)
      break;
  }
}

/***********************************************************************
 * Run beta-binoimial regression using the specified table with
 * proportions and design matrix
 */
int
main_radmeth(int argc, char *argv[]) {
  try {
    static const std::string description =
      "calculate differential methylation scores";

    std::string outfile;
    std::string test_factor;
    bool verbose{false};
    bool show_progress{false};
    bool more_na_info{false};
    std::uint32_t n_threads{1};

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), description,
                           "<design-matrix> <data-matrix>");
    opt_parse.set_show_defaults();
    opt_parse.add_opt("out", 'o', "output file", true, outfile);
    opt_parse.add_opt("threads", 't', "number of threads to use", false,
                      n_threads);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, verbose);
    opt_parse.add_opt("progress", '\0', "show progress", false, show_progress);
    opt_parse.add_opt(
      "na-info", 'n',
      "if a p-value is not calculated, print NAs in more detail: "
      "low count (NA_LOW_COV) extreme values (NA_EXTREME_CNT) or "
      "numerical errors in likelihood ratios (NA)",
      false, more_na_info);
    opt_parse.add_opt("factor", 'f', "a factor to test", true, test_factor);

    std::vector<std::string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      std::cerr << opt_parse.help_message() << std::endl
                << opt_parse.about_message() << std::endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      std::cerr << opt_parse.about_message() << std::endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      std::cerr << opt_parse.option_missing_message() << std::endl;
      return EXIT_SUCCESS;
    }
    if (leftover_args.size() != 2) {
      std::cerr << opt_parse.help_message() << std::endl;
      return EXIT_SUCCESS;
    }
    const std::string design_filename(leftover_args.front());
    const std::string table_filename(leftover_args.back());
    /****************** END COMMAND LINE OPTIONS *****************/

    // initialize full design matrix from file
    Regression orig_alt_model;
    read_design(verbose, design_filename, orig_alt_model.design);

    // Check that provided test factor name exists and find its index.
    const auto test_factor_idx =
      get_test_factor_idx(orig_alt_model, test_factor);

    // verify that the design includes more than one level for the
    // test factor
    if (!has_two_values(orig_alt_model, test_factor_idx)) {
      const auto first_level = orig_alt_model.design.matrix[0][test_factor_idx];
      throw std::runtime_error("test factor only one level: " + test_factor +
                               ", " + std::to_string(first_level));
    }

    Regression orig_null_model;
    orig_null_model.design = orig_alt_model.design;
    remove_factor(orig_null_model.design, test_factor_idx);

    radmeth(show_progress, more_na_info, n_threads, table_filename, outfile,
            orig_alt_model, orig_null_model, test_factor_idx);
  }
  catch (const std::exception &e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
