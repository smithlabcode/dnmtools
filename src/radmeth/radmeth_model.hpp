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

#include <cstdint>
#include <istream>
#include <ostream>
#include <string>
#include <vector>

struct mcounts {
  std::uint32_t n_reads{};
  std::uint32_t n_meth{};
};

[[nodiscard]] inline std::istream &
operator>>(std::istream &in, mcounts &rm) {
  return in >> rm.n_reads >> rm.n_meth;
}

struct cumul_counts {
  std::vector<std::uint32_t> m_counts;
  std::vector<std::uint32_t> r_counts;
  std::vector<std::uint32_t> d_counts;
};

struct Regression {
  static double tolerance;        // 1e-3;
  static double stepsize;         // 0.001;
  static std::uint32_t max_iter;  // 250;

  Design design;
  std::string rowname;
  std::vector<mcounts> mc;
  double max_loglik{};

  // scratch space
  std::vector<cumul_counts> cumul;
  std::vector<double> p_v;
  std::vector<double> cache;
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

#endif  // RADMETH_MODEL_HPP
