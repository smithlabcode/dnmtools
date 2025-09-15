/* Copyright (C) 2013-2025 Andrew D Smith
 *
 * Author: Andrew D. Smith
 * Contributors: Egor Dolzhenko and Guilherme Sena
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

#include "radmeth_model.hpp"
#include "radmeth_optimize.hpp"

// smithlab headers
#include "GenomicRegion.hpp"
#include "OptionParser.hpp"
#include "smithlab_os.hpp"
#include "smithlab_utils.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

[[nodiscard]] static std::string
format_duration(const std::chrono::duration<double> elapsed) {
  static constexpr auto s_per_h = 3600;
  static constexpr auto s_per_m = 60;
  const double tot_s = elapsed.count();

  // break down into hours, minutes, seconds
  const std::uint32_t hours = tot_s / 3600;
  const std::uint32_t minutes = (static_cast<int>(tot_s) % s_per_h) / s_per_m;
  const double seconds = tot_s - (hours * s_per_h) - (minutes * s_per_m);

  std::ostringstream oss;
  oss << std::setfill('0') << std::setw(2) << hours << ":" << std::setfill('0')
      << std::setw(2) << minutes << ":" << std::fixed << std::setprecision(2)
      << std::setw(5) << seconds;
  return oss.str();
}

struct file_progress {
  double one_thousand_over_filesize{};
  std::size_t prev_offset{};
  explicit file_progress(const std::string &filename) :
    one_thousand_over_filesize{1000.0 / std::filesystem::file_size(filename)} {}
  void
  operator()(std::ifstream &in) {
    const std::size_t curr_offset =
      in.eof() ? 1000 : in.tellg() * one_thousand_over_filesize;
    if (curr_offset <= prev_offset)
      return;
    std::ios old_state(nullptr);
    old_state.copyfmt(std::cerr);
    std::cerr << "\r[progress: " << std::setw(5) << std::fixed
              << std::setprecision(1) << (curr_offset / 10.0)
              << (curr_offset == 1000 ? "%]\n" : "%]");
    std::cerr.copyfmt(old_state);
    prev_offset = (curr_offset == 1000)
                    ? std::numeric_limits<std::size_t>::max()
                    : curr_offset;
  }
};

[[nodiscard]] static bool
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

[[nodiscard]] static std::vector<std::string>
get_sample_names_from_header(const std::string &header) {
  std::istringstream iss(header);
  std::string token;
  std::vector<std::string> sample_names;
  while (iss >> token)
    sample_names.push_back(token);
  return sample_names;
}

// Series representation for the lower incomplete gamma P(a,x)
[[nodiscard]] static double
gamma_p_series(const double a, const double x) {
  constexpr auto eps = std::numeric_limits<double>::epsilon();
  constexpr auto max_iter = 100;
  double sum = 1.0 / a;
  double term = sum;
  for (auto n = 1; n < max_iter; ++n) {
    term *= x / (a + n);
    sum += term;
    if (term < eps * sum)
      break;
  }
  return sum * std::exp(-x + a * std::log(x) - std::lgamma(a));
}

[[nodiscard]] static inline double
safe_floor(const double x, const double floor_val) {
  return std::abs(x) < floor_val ? floor_val : x;
}

// Continued fraction representation for the upper incomplete gamma Q(a,x)
[[nodiscard]] static double
gamma_q_contfrac(const double a, const double x) {
  constexpr auto epsilon = std::numeric_limits<double>::epsilon();
  constexpr auto fpmin = std::numeric_limits<double>::min() / epsilon;
  constexpr auto max_iter = 100;
  double b = x + 1 - a;
  double c = 1.0 / fpmin;
  double d = 1.0 / b;
  double h = d;
  for (auto i = 1; i < max_iter; ++i) {
    const double an = -i * (i - a);
    b += 2;
    d = safe_floor(an * d + b, fpmin);
    c = safe_floor(b + an / c, fpmin);
    d = 1.0 / d;
    const double delta = d * c;
    h *= delta;
    if (std::abs(delta - 1.0) < epsilon)
      break;
  }
  return std::exp(-x + a * std::log(x) - std::lgamma(a)) * h;
}

// Regularized lower incomplete gamma P(a,x)
[[nodiscard]] static double
gamma_p(const double a, const double x) {
  if (x < 0 || a <= 0)
    return 0.0;
  if (x == 0)
    return 0.0;
  if (x < a + 1.0)
    return gamma_p_series(a, x);
  return 1.0 - gamma_q_contfrac(a, x);
}

// chi-square CDF: P(k/2, x/2)
[[nodiscard]] static double
chi_square_cdf(const double x, const double k) {
  return gamma_p(k * 0.5, x * 0.5);
}

// Given the maximum likelihood estimates of the full and reduced models, the
// function outputs the p-value of the log-likelihood ratio. *Note* that it is
// assumed that the reduced model has one fewer factor than the reduced model.
[[nodiscard]] static double
llr_test(const double null_loglik, const double full_loglik) {
  // The log-likelihood ratio statistic.
  const double log_lik_stat = -2 * (null_loglik - full_loglik);

  // It is assumed that null model has one fewer factor than the full
  // model. Hence the number of degrees of freedom is 1.
  const std::size_t degrees_of_freedom = 1;

  // Log-likelihood ratio statistic has a chi-sqare distribution.
  const double chisq_p = chi_square_cdf(log_lik_stat, degrees_of_freedom);

  const double p_value = 1.0 - chisq_p;

  return p_value;
}

static bool
has_low_coverage(const Regression &reg, const std::size_t test_factor) {
  bool cvrd_in_test_fact_smpls = false;
  const auto &tcol = reg.design.tmatrix[test_factor];
  for (std::size_t i = 0; i < reg.n_samples() && !cvrd_in_test_fact_smpls; ++i)
    cvrd_in_test_fact_smpls = (tcol[i] == 1 && reg.props.mc[i].n_reads != 0);

  bool cvrd_in_other_smpls = false;
  for (std::size_t i = 0; i < reg.n_samples() && !cvrd_in_other_smpls; ++i)
    cvrd_in_other_smpls = (tcol[i] != 1 && reg.props.mc[i].n_reads != 0);

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
  const auto &tcol = reg.design.tmatrix[test_factor];
  for (const auto x : tcol)
    if (x != tcol[0])
      return true;
  return false;
}

[[nodiscard]] static std::uint32_t
get_test_factor_idx(const Regression &model, const std::string &test_factor) {
  const auto &factors = model.design.factor_names;
  const auto itr =
    std::find(std::cbegin(factors), std::cend(factors), test_factor);

  if (itr == std::cend(factors))
    throw std::runtime_error("factor not part of design: " + test_factor);

  return std::distance(std::cbegin(factors), itr);
}

[[nodiscard]] static Design
read_design(const std::string &design_filename) {
  std::ifstream design_file(design_filename);
  if (!design_file)
    throw std::runtime_error("could not open file: " + design_filename);
  Design design;
  design_file >> design;
  return design;
}

enum class row_status : std::uint8_t {
  ok,
  na,
  na_low_cov,
  na_extreme_cnt,
};

static void
radmeth(const bool show_progress, const bool more_na_info,
        const std::uint32_t n_threads, const std::string &table_filename,
        const std::string &outfile, const Regression &alt_model,
        const Regression &null_model, const std::uint32_t test_factor_idx) {
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
  if (!consistent_sample_names(alt_model, sample_names_header)) {
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

  std::vector<Regression> alt_models(n_threads, alt_model);
  std::vector<Regression> null_models(n_threads, null_model);

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
      threads.emplace_back([&, thread_id] {
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
          t_alt_model.props.parse(lines[b]);
          if (t_alt_model.props_size() != n_samples)
            throw std::runtime_error("found row with wrong number of columns");

          const auto [p_val, status] = [&]() -> std::tuple<double, row_status> {
            // Skip the test if (1) no coverage in all cases or in all
            // controls, or (2) the site is completely methylated or
            // completely unmethylated across all samples.
            if (has_low_coverage(t_alt_model, test_factor_idx))
              return std::tuple{1.0, row_status::na_low_cov};

            if (has_extreme_counts(t_alt_model))
              return std::tuple{1.0, row_status::na_extreme_cnt};

            fit_regression_model(t_alt_model, p_estim_alt, phi_estim_alt);

            t_null_model.props = t_alt_model.props;

            fit_regression_model(t_null_model, p_estim_null, phi_estim_null);

            const double p_value =
              llr_test(t_null_model.max_loglik, t_alt_model.max_loglik);

            return (p_value != p_value) ? std::tuple{1.0, row_status::na}
                                        : std::tuple{p_value, row_status::ok};
          }();

          n_bytes[b] = [&] {
            // clang-format off
            const int n_prefix_bytes =
              std::sprintf(bufs[b].data(), "%s\t",
                           t_alt_model.props.rowname.data());
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

            const int n_param_bytes = [&] {
              std::int32_t n_bytes = 0;
              if (p_estim_alt.empty())
                p_estim_alt.resize(n_groups, 0.0);
              for (auto g_idx = 0u; g_idx < n_groups; ++g_idx) {
                const int n = std::sprintf(cursor, "\t%f", p_estim_alt[g_idx]);
                cursor += n;
                if (n < 0)
                  return n;
                n_bytes += n;
              }
              const int n = std::sprintf(cursor, "\t%f\n", phi_estim_alt);
              if (n < 0)
                return n;
              n_bytes += n;
              return n_bytes;
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

static void
ensure_sample_order(const std::string &table_filename, Regression &alt_model,
                    Regression &null_model) {
  std::ifstream table_file(table_filename);
  if (!table_file)
    throw std::runtime_error("could not open file: " + table_filename);
  std::string header;
  std::getline(table_file, header);
  const auto sample_names = get_sample_names_from_header(header);
  alt_model.order_samples(sample_names);
  null_model.order_samples(sample_names);
}

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
    opt_parse.add_opt("tolerance", '\0',
                      "numerical tolerance to test convergence", false,
                      Regression::tolerance);
    opt_parse.add_opt("max-iter", '\0',
                      "max iterations when estimating parameters", false,
                      Regression::max_iter);

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

    if (verbose)
      std::cerr << "design table filename: " << design_filename << "\n\n";

    // initialize full design matrix from file
    Regression alt_model;
    alt_model.design = read_design(design_filename);

    if (verbose)
      std::cerr << "Alternate model:\n" << alt_model.design << '\n';

    // Check that provided test factor name exists and find its index.
    const auto test_factor_idx = get_test_factor_idx(alt_model, test_factor);

    // verify that the design includes more than one level for the
    // test factor
    if (!has_two_values(alt_model, test_factor_idx)) {
      const auto first_level = alt_model.design.matrix[0][test_factor_idx];
      throw std::runtime_error("test factor only one level: " + test_factor +
                               ", " + std::to_string(first_level));
    }

    Regression null_model = alt_model;
    null_model.design = alt_model.design.drop_factor(test_factor_idx);
    if (verbose)
      std::cerr << "Null model:\n" << null_model.design << '\n';

    ensure_sample_order(table_filename, alt_model, null_model);

    // clang-format off
    if (verbose)
      std::cerr << "Output columns:\n"
                << "(1) chrom\n"
                << "(2) position\n"
                << "(3) strand\n"
                << "(4) cytosine context\n"
                << "(5) p-value\n"
                << "(6) estimated methylation 0\n"
                << "(7) estimated methylation 1\n"
                << "(8) estimated dispersion\n"
                << "estimated methylation is for test factor value (0 or 1)\n"
                << '\n';
    // clang-format on
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
