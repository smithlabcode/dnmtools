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

#ifndef RADMETH_MODEL_HPP
#define RADMETH_MODEL_HPP

#include <cstdint>
#include <istream>
#include <ostream>
#include <string>
#include <vector>

struct Design {
  std::vector<std::string> factor_names;
  std::vector<std::string> sample_names;
  std::vector<std::vector<std::uint8_t>> matrix;
  std::vector<std::vector<std::uint8_t>> tmatrix;
  std::size_t
  n_factors() const {
    return factor_names.size();
  }
  std::size_t
  n_samples() const {
    return sample_names.size();
  }
};

std::istream &
operator>>(std::istream &is, Design &design);

std::ostream &
operator<<(std::ostream &os, const Design &design);

struct mcounts {
  std::uint32_t n_reads{};
  std::uint32_t n_meth{};
};

[[nodiscard]] inline std::istream &
operator>>(std::istream &in, mcounts &rm) {
  return in >> rm.n_reads >> rm.n_meth;
}

struct SiteProportions {
  std::string chrom;
  std::size_t position{};
  char strand{};
  std::string context;
  std::vector<mcounts> mc;

  void
  parse(const std::string &line);
};

struct Regression {
  static double tolerance;        // 1e-4;
  static double stepsize;         // 0.001;
  static std::uint32_t max_iter;  // 700;

  Design design;
  SiteProportions props;
  double max_loglik{};

  [[nodiscard]] std::size_t
  n_factors() const {
    return design.n_factors();
  }

  [[nodiscard]] std::size_t
  props_size() const {
    return std::size(props.mc);
  }

  [[nodiscard]] std::size_t
  n_samples() const {
    return design.n_samples();
  }
};

#endif
