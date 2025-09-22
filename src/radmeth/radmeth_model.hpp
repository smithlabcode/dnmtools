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

#include <cassert>
#include <cstdint>
#include <istream>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

struct cumul_counts {
  std::vector<std::uint32_t> m_counts;
  std::vector<std::uint32_t> r_counts;
  std::vector<std::uint32_t> d_counts;
};

struct Design {
  std::vector<std::string> factor_names;
  std::vector<std::string> sample_names;
  std::vector<std::vector<std::uint8_t>> matrix;   // samples=rows, factors=cols
  std::vector<std::vector<std::uint8_t>> tmatrix;  // factors=rows, samples=cols
  std::vector<std::vector<std::uint8_t>> groups;   // combs of fact levels
  std::vector<std::uint32_t> group_id;             // map sample to group
  std::vector<std::uint32_t> g_idx_map;  // larger group id into groups

  [[nodiscard]] std::uint32_t
  n_factors() const {
    return std::size(factor_names);
  }

  [[nodiscard]] std::uint32_t
  n_groups() const {
    return std::size(groups);
  }

  [[nodiscard]] std::uint32_t
  n_params() const {
    return std::size(groups);
  }

  [[nodiscard]] std::uint32_t
  n_samples() const {
    return std::size(sample_names);
  }

  [[nodiscard]] std::uint32_t
  get_g_idx(const std::uint32_t g_idx) const {
    if (g_idx >= std::size(g_idx_map))
      throw std::runtime_error("bad g_idx value: " + std::to_string(g_idx));
    return g_idx_map[g_idx];
  }

  [[nodiscard]] Design
  drop_factor(const std::uint32_t factor_idx);

  void
  order_samples(const std::vector<std::string> &ordered_names);
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

struct SiteProp {
  std::string rowname;
  std::vector<mcounts> mc;

  void
  parse(const std::string &line);
};

struct Regression {
  static double tolerance;        // 1e-3;
  static double stepsize;         // 0.001;
  static std::uint32_t max_iter;  // 250;

  Design design_p;
  Design design_disp;
  SiteProp props;
  double max_loglik{};

  // scratch space
  std::vector<cumul_counts> cumul;
  std::vector<double> p_v;
  std::vector<std::vector<double>> disp_cache;
  std::vector<std::uint32_t> max_r_count;  // {};

  [[nodiscard]] std::uint32_t
  n_disp_factors() const {
    return design_disp.n_factors();
  }

  [[nodiscard]] std::uint32_t
  n_disp_groups() const {
    return design_disp.n_groups();
  }

  [[nodiscard]] std::uint32_t
  n_disp_params() const {
    return design_disp.n_params();
  }

  [[nodiscard]] std::uint32_t
  get_disp_idx(const std::uint32_t g_idx) const {
    return design_disp.get_g_idx(g_idx);
  }

  [[nodiscard]] std::uint32_t
  get_disp_param_idx(const std::uint32_t idx) const {
    return design_p.n_params() + idx;
  }

  [[nodiscard]] const std::vector<std::uint8_t> &
  get_disp_group(const std::uint32_t g_idx) const {
    return design_disp.groups[get_disp_idx(g_idx)];
  }

  [[nodiscard]] std::uint32_t
  n_p_factors() const {
    return design_disp.n_factors();
  }

  [[nodiscard]] std::uint32_t
  n_p_groups() const {
    return design_disp.n_groups();
  }

  [[nodiscard]] std::uint32_t
  n_p_params() const {
    return design_disp.n_params();
  }

  [[nodiscard]] std::uint32_t
  get_p_idx(const std::uint32_t g_idx) const {
    return design_p.get_g_idx(g_idx);
  }

  [[nodiscard]] const std::vector<std::uint8_t> &
  get_p_group(const std::uint32_t g_idx) const {
    return design_p.groups[get_p_idx(g_idx)];
  }

  [[nodiscard]] std::uint32_t
  props_size() const {
    return std::size(props.mc);
  }

  [[nodiscard]] std::uint32_t
  n_samples() const {
    assert(design_p.n_samples() == design_disp.n_samples());
    return design_p.n_samples();
  }

  [[nodiscard]] std::uint32_t
  n_groups() const {
    return std::max(design_p.n_groups(), design_disp.n_groups());
  }

  [[nodiscard]] std::uint32_t
  n_params() const {
    return n_p_params() + n_disp_params();
  }

  [[nodiscard]] const std::vector<std::uint32_t> &
  get_group_id() const {
    return design_p.n_groups() > design_disp.n_groups() ? design_p.group_id
                                                        : design_disp.group_id;
  }

  [[nodiscard]] const std::vector<std::string> &
  get_sample_names() const {
    assert(design_p.sample_names == design_disp.sample_names);
    return design_p.sample_names;
  }

  [[nodiscard]] const std::vector<std::vector<std::uint8_t>> &
  get_tmatrix() const {
    return design_p.n_groups() > design_disp.n_groups() ? design_p.tmatrix
                                                        : design_disp.tmatrix;
  }

  void
  order_samples(const std::vector<std::string> &ordered_names) {
    design_disp.order_samples(ordered_names);
    design_p.order_samples(ordered_names);
  }
};

#endif
