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

#ifndef RADMETH_DESIGN_HPP
#define RADMETH_DESIGN_HPP

#include <cstdint>
#include <istream>
#include <ostream>
#include <string>
#include <vector>

struct Design {
  std::vector<std::string> factor_names;
  std::vector<std::string> sample_names;
  std::vector<std::vector<std::uint8_t>> matrix;   // samples=rows, factors=cols
  std::vector<std::vector<std::uint8_t>> tmatrix;  // factors=rows, samples=cols
  std::vector<std::vector<std::uint8_t>> groups;   // combs of fact levels
  std::vector<std::uint32_t> group_id;             // assign group to sample

  [[nodiscard]] static Design
  read_design(const std::string &design_filename);

  [[nodiscard]] std::size_t
  n_factors() const {
    return std::size(factor_names);
  }

  [[nodiscard]] std::size_t
  n_groups() const {
    return std::size(groups);
  }

  [[nodiscard]] std::size_t
  n_samples() const {
    return std::size(sample_names);
  }

  [[nodiscard]] Design
  drop_factor(const std::size_t factor_idx);

  void
  order_samples(const std::vector<std::string> &ordered_names);

  [[nodiscard]] std::uint32_t
  get_test_factor_idx(const std::string &test_factor) const;

  [[nodiscard]] bool
  has_two_values(const std::size_t test_factor) const;
};

std::istream &
operator>>(std::istream &is, Design &design);

std::ostream &
operator<<(std::ostream &os, const Design &design);

void
ensure_sample_order(const std::string &table_filename, Design &design);

[[nodiscard]] std::vector<std::string>
get_sample_names_from_header(const std::string &header);

[[nodiscard]] bool
consistent_sample_names(const Design &design, const std::string &header);

#endif  // RADMETH_DESIGN_HPP
