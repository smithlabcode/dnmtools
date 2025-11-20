/* Copyright (C) 2025 Andrew D Smith
 *
 * Author: Andrew D Smith
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

#ifndef RADMETH_MODEL_HPP
#define RADMETH_MODEL_HPP

#include "radmeth_design.hpp"

#include <algorithm>
#include <charconv>
#include <cstdint>
#include <istream>
#include <ostream>
#include <string>
#include <type_traits>
#include <vector>

template <typename T> struct mcounts {
  std::uint32_t n_reads{};
  T n_meth{};
};

template <typename T>
[[nodiscard]] inline std::istream &
operator>>(std::istream &in, mcounts<T> &rm) {
  return in >> rm.n_reads >> rm.n_meth;
}

struct cumul_counts {
  std::vector<std::uint32_t> m_counts;
  std::vector<std::uint32_t> r_counts;
  std::vector<std::uint32_t> d_counts;
};

struct vars_cache {
  double p{};
  double a{};
  double b{};
  double lgamma_a{};
  double lgamma_b{};
  double lgamma_a_b{};
  double digamma_a{};
  double digamma_b{};
  double digamma_a_b{};
};

template <typename T> struct Regression {
  static double tolerance;        // 1e-3;
  static double stepsize;         // 0.001;
  static std::uint32_t max_iter;  // 250;
  static constexpr double pseudocount = 1.0 / 256.0;

  Design design;
  std::string rowname;
  std::vector<mcounts<T>> mc;
  double max_loglik{};

  // scratch space
  std::vector<cumul_counts> cumul;
  std::vector<double> p_v;
  std::vector<double> cache_log1p_factors;
  std::vector<double> cache_dispersion_effect;
  std::vector<vars_cache> cache;  // scratch space
  std::uint32_t max_r_count{};

  void
  parse(const std::string &line);

  [[nodiscard]] std::size_t
  n_factors() const {
    return design.n_factors();
  }

  [[nodiscard]] std::size_t
  n_groups() const {
    return design.n_groups();
  }

  [[nodiscard]] std::size_t
  n_params() const {
    return n_factors() + 1;
  }

  [[nodiscard]] std::size_t
  props_size() const {
    return std::size(mc);
  }

  [[nodiscard]] std::size_t
  n_samples() const {
    return design.n_samples();
  }
};

template <typename T>
void
Regression<T>::parse(const std::string &line) {
  const auto first_ws = line.find_first_of(" \t");

  // Parse the row name (must be like: "chr:position:strand:context")
  bool failed = (first_ws == std::string::npos);

  auto field_s = line.data();
  auto field_e = line.data() + first_ws;
  rowname = std::string{field_s, field_e};
  std::replace(std::begin(rowname), std::end(rowname), ':', '\t');
  if (failed)
    throw std::runtime_error("failed to parse label from:\n" + line);

  // Parse the counts of total reads and methylated reads
  const auto is_sep = [](const auto x) { return std::isspace(x); };
  const auto not_sep = [](const auto x) { return std::isdigit(x); };

  const auto line_end = line.data() + std::size(line);
  mcounts<T> mc1{};
  mc.clear();
  while (field_e != line_end) {
    // get the total count
    field_s = std::find_if(field_e + 1, line_end, not_sep);
    field_e = std::find_if(field_s + 1, line_end, is_sep);
    {
      const auto [ptr, ec] = std::from_chars(field_s, field_e, mc1.n_reads);
      failed = failed || (ec != std::errc{});
    }

    field_s = std::find_if(field_e + 1, line_end, not_sep);
    field_e = std::find_if(field_s + 1, line_end, is_sep);

    if constexpr (std::is_floating_point_v<T>) {
#ifdef __APPLE__
      const int ret = std::sscanf(field_s, "%lf", &mc1.n_meth);
      failed = failed || (ret < 1);
#else
      const auto [ptr, ec] = std::from_chars(field_s, field_e, mc1.n_meth);
      failed = failed || (ec != std::errc{});
#endif
    }
    else {
      // Apple clang can do std::from_chars for int types
      const auto [ptr, ec] = std::from_chars(field_s, field_e, mc1.n_meth);
      failed = failed || (ec != std::errc{});
    }

    mc.push_back(mc1);
  }

  const auto add_pseudocount = [&](auto x) {
    x.n_meth = x.n_reads > 0.0 ? pseudocount / 2.0 +
                                   x.n_meth * (1.0 - pseudocount / x.n_reads)
                               : 0.0;
    return x;
  };

  std::transform(std::begin(mc), std::end(mc), std::begin(mc), add_pseudocount);

  if (failed)
    throw std::runtime_error("failed to parse counts from:\n" + line);
}

#endif  // RADMETH_MODEL_HPP
