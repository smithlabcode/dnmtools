/* Copyright (C) 2025 Andrew D Smith
 *
 * Author: Andrew D. Smith
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

#include "radmeth_utils.hpp"

#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <sstream>
#include <string>

[[nodiscard]] std::string
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

file_progress::file_progress(const std::string &filename) :
  one_thousand_over_filesize{1000.0 / std::filesystem::file_size(filename)} {}

void
file_progress::operator()(std::ifstream &in) {
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
  prev_offset = (curr_offset == 1000) ? std::numeric_limits<std::size_t>::max()
                                      : curr_offset;
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
[[nodiscard]] double
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
