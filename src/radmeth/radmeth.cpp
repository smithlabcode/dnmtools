/* Copyright (C) 2013-2025 Andrew D Smith, Egor Dolzhenko and Guilherme Sena
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
 */

#include "file_progress.hpp"
#include "radmeth_design.hpp"
#include "radmeth_model.hpp"
#include "radmeth_optimize_params.hpp"
#include "radmeth_optimize_series.hpp"
#include "radmeth_utils.hpp"

#include "OptionParser.hpp"

#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <thread>
#include <tuple>
#include <utility>
#include <vector>

template <typename RegressionType>
[[nodiscard]] static bool
has_low_coverage(const RegressionType &reg, const std::size_t test_factor) {
  bool cvrd_in_test_fact_smpls = false;
  const auto &tcol = reg.design.tmatrix[test_factor];
  for (std::size_t i = 0; i < reg.n_samples() && !cvrd_in_test_fact_smpls; ++i)
    cvrd_in_test_fact_smpls = (tcol[i] == 1 && reg.mc[i].n_reads != 0);

  bool cvrd_in_other_smpls = false;
  for (std::size_t i = 0; i < reg.n_samples() && !cvrd_in_other_smpls; ++i)
    cvrd_in_other_smpls = (tcol[i] != 1 && reg.mc[i].n_reads != 0);

  return !cvrd_in_test_fact_smpls || !cvrd_in_other_smpls;
}

template <typename RegressionType>
[[nodiscard]] static bool
has_extreme_counts(const RegressionType &reg) {
  const auto &mc = reg.mc;

  bool full_meth = true;
  for (auto i = std::cbegin(mc); i != std::cend(mc) && full_meth; ++i)
    full_meth = (i->n_reads == i->n_meth);

  bool no_meth = true;
  for (auto i = std::cbegin(mc); i != std::cend(mc) && no_meth; ++i)
    no_meth = (i->n_meth == 0.0);

  return full_meth || no_meth;
}

enum class row_status : std::uint8_t {
  ok,
  na,
  na_low_cov,
  na_extreme_cnt,
};

template <typename T>
static void
radmeth(const bool show_progress, const bool more_na_info,
        const std::uint32_t n_threads, const std::string &table_filename,
        const std::string &outfile, const Regression<T> &alt_model,
        const Regression<T> &null_model, const std::uint32_t test_factor_idx) {
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

  const auto sample_names = get_sample_names_from_header(sample_names_header);
  if (!consistent_sample_names(alt_model.design, sample_names_header)) {
    // clang-format off
    const auto message =
      R"(
does not match factor names or their order in the design matrix. Check
that the design matrix and the proportion table are correctly formatted.
)";
    // clang-format on
    throw std::runtime_error("header:\n" + sample_names_header + message);
  }

  const std::size_t n_samples = alt_model.n_samples();

  std::ofstream out(outfile);
  if (!out)
    throw std::runtime_error("failed to open output file: " + outfile);

  file_progress progress{table_filename};

  std::vector<std::vector<char>> bufs(max_lines,
                                      std::vector<char>(buf_size, 0));
  std::vector<int> n_bytes(max_lines, 0);
  std::vector<std::string> lines(max_lines);

  std::vector<Regression<T>> alt_models(n_threads, alt_model);
  std::vector<Regression<T>> null_models(n_threads, null_model);

  const auto n_groups = alt_model.n_groups();

  // Iterate over rows in the file and do llr test on proportions from
  // each. Do this in sets of rows to avoid having to spawn too many threads.
  while (true) {
    std::uint32_t n_lines = 0;
    while (n_lines < max_lines && std::getline(table_file, lines[n_lines]))
      ++n_lines;
    if (n_lines == 0)
      break;

    std::vector<std::thread> threads;
    for (auto thread_id = 0u; thread_id < n_threads; ++thread_id) {
      threads.emplace_back(  // NOLINT(performance-inefficient-vector-operation)
        [&, thread_id] {
          std::vector<double> p_estim_alt;
          std::vector<double> p_estim_null;
          double phi_estim_alt{};
          double phi_estim_null{};

          auto &t_alt_model = alt_models[thread_id];
          auto &t_null_model = null_models[thread_id];
          for (auto b = 0u; b < n_lines; ++b) {
            // ADS: rows done by different threads are interleaved because the
            // difficult (e.g., high-coverage) rows can be consecutive and this
            // balances work better.
            if (b % n_threads != thread_id)
              continue;
            t_alt_model.parse(lines[b]);
            if (t_alt_model.props_size() != n_samples)
              throw std::runtime_error(
                "found row with wrong number of columns");

            const auto p_val_status = [&]() -> std::tuple<double, row_status> {
              // Skip the test if (1) no coverage in all cases or in all
              // controls, or (2) the site is completely methylated or
              // completely unmethylated across all samples.
              if (has_low_coverage(t_alt_model, test_factor_idx))
                return std::tuple{1.0, row_status::na_low_cov};

              if (has_extreme_counts(t_alt_model))
                return std::tuple{1.0, row_status::na_extreme_cnt};

              fit_regression_model(t_alt_model, p_estim_alt, phi_estim_alt);

              t_null_model.mc = t_alt_model.mc;
              t_null_model.rowname = t_alt_model.rowname;

              fit_regression_model(t_null_model, p_estim_null, phi_estim_null);

              const double p_value =
                llr_test(t_null_model.max_loglik, t_alt_model.max_loglik);

              return (p_value != p_value) ? std::tuple{1.0, row_status::na}
                                          : std::tuple{p_value, row_status::ok};
            }();
            // ADS: avoid capture structured binding in C++17
            const auto p_val = std::get<0>(p_val_status);
            const auto status = std::get<1>(p_val_status);

            n_bytes[b] = [&] {
              auto bufsize = std::size(bufs[b]);
              // clang-format off
            const int n_prefix_bytes =
              std::snprintf(bufs[b].data(), bufsize,
                           "%s\t", t_alt_model.rowname.data());
              // clang-format on
              if (n_prefix_bytes < 0)
                return n_prefix_bytes;

              bufsize -= n_prefix_bytes;
              // NOLINTNEXTLINE(*-pointer-arithmetic)
              auto cursor = bufs[b].data() + n_prefix_bytes;

              const int n_pval_bytes = [&] {
                if (status == row_status::ok)
                  return std::snprintf(cursor, bufsize, "%.6g", p_val);
                if (!more_na_info || status == row_status::na)
                  return std::snprintf(cursor, bufsize, "NA");
                if (status == row_status::na_extreme_cnt)
                  return std::snprintf(cursor, bufsize, "NA_EXTREME_CNT");
                // if (status == row_status::na_low_cov)
                return std::snprintf(cursor, bufsize, "NA_LOW_COV");
              }();

              if (n_pval_bytes < 0)
                return n_pval_bytes;

              bufsize -= n_pval_bytes;
              cursor += n_pval_bytes;  // NOLINT(*-pointer-arithmetic)

              const int n_param_bytes = [&] {
                std::int32_t n_param_bytes = 0;
                if (p_estim_alt.empty())
                  p_estim_alt.resize(n_groups, 0.0);
                for (auto g_idx = 0u; g_idx < n_groups; ++g_idx) {
                  const int n =
                    std::snprintf(cursor, bufsize, "\t%f", p_estim_alt[g_idx]);
                  bufsize -= n;
                  cursor += n;  // NOLINT(*-pointer-arithmetic)
                  if (n < 0)
                    return n;
                  n_param_bytes += n;
                }
                const auto od = overdispersion_factor(n_samples, phi_estim_alt);
                const int n = std::snprintf(cursor, bufsize, "\t%f\n", od);
                if (n < 0)
                  return n;
                n_param_bytes += n;
                return n_param_bytes;
              }();

              if (n_param_bytes < 0)
                return n_param_bytes;

              return n_prefix_bytes + n_pval_bytes + n_param_bytes;
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

  if (show_progress)
    progress(table_file);
}

int
main_radmeth(int argc, char *argv[]) {  // NOLINT(*-avoid-c-arrays)
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
    OptionParser opt_parse(argv[0],  // NOLINT(*-pointer-arithmetic)
                           description, "<design-matrix> <data-matrix>");
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
    opt_parse.add_opt("tolerance", '\0',
                      "numerical tolerance to test convergence", false,
                      radmeth_optimize_params::tolerance);
    opt_parse.add_opt("max-iter", '\0',
                      "max iterations when estimating parameters", false,
                      radmeth_optimize_params::max_iter);
    std::vector<std::string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      std::cerr << opt_parse.help_message() << '\n'
                << opt_parse.about_message() << '\n';
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
    if (leftover_args.size() != 2) {
      std::cerr << opt_parse.help_message() << '\n';
      return EXIT_SUCCESS;
    }
    const std::string design_filename(leftover_args.front());
    const std::string table_filename(leftover_args.back());
    /****************** END COMMAND LINE OPTIONS *****************/

    Design design = Design::read_design(design_filename);

    // Check that provided test factor name exists and find its index.
    const auto test_factor_idx = design.get_test_factor_idx(test_factor);

    ensure_sample_order(table_filename, design);

    // verify that the design includes more than one level for the
    // test factor
    if (!design.has_two_values(test_factor_idx)) {
      const auto first_level = design.matrix[0][test_factor_idx];
      throw std::runtime_error("test factor only one level: " + test_factor +
                               ", " + std::to_string(first_level));
    }

    const Design null_design = design.drop_factor(test_factor_idx);

    // clang-format off
    if (verbose)
      std::cerr << "design table filename: " << design_filename << "\n\n"
                << "Alternate model:\n" << design << '\n'
                << "Null model:\n" << null_design << '\n'
                << "Output columns:\n"
                << "(1) chrom\n"
                << "(2) position\n"
                << "(3) strand\n"
                << "(4) cytosine context\n"
                << "(5) p-value\n"
                << "(6) overdispersion factor (variance above binomial)\n"
                << "(7-) estimated methylation for each unique group\n"
                << "groups are defined by combinations of factor levels\n\n";
    // clang-format on

    Regression<std::uint32_t> alt_model;
    alt_model.design = design;
    Regression<std::uint32_t> null_model;
    null_model.design = null_design;

    const auto start_time = std::chrono::steady_clock::now();
    radmeth(show_progress, more_na_info, n_threads, table_filename, outfile,
            alt_model, null_model, test_factor_idx);
    const auto stop_time = std::chrono::steady_clock::now();

    if (verbose)
      std::cerr << "[total time: " << format_duration(stop_time - start_time)
                << "]\n";
  }
  catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
