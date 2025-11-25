/* Copyright (C) 2025 Andrew D Smith
 *
 * Authors: Andrew D. Smith
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

#include "radmeth_design.hpp"

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <istream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

static void
make_groups(Design &design) {
  auto s = design.matrix;
  std::sort(std::begin(s), std::end(s));
  s.erase(std::unique(std::begin(s), std::end(s)), std::end(s));
  design.groups = std::move(s);
}

static void
assign_groups(Design &design) {
  const auto &matrix = design.matrix;
  const auto &s = design.groups;
  const auto n_samples = design.n_samples();
  auto &group_id = design.group_id;
  group_id.resize(n_samples);
  for (auto i = 0u; i < n_samples; ++i) {
    const auto x = std::find(std::cbegin(s), std::cend(s), matrix[i]);
    group_id[i] = std::distance(std::cbegin(s), x);
  }
}

template <typename T>
static void
transpose(const std::vector<std::vector<T>> &mat,
          std::vector<std::vector<T>> &tmat) {
  const auto n_row = std::size(mat);
  const auto n_col = std::size(mat.front());
  tmat.resize(n_col, std::vector<T>(n_row, 0.0));
  for (auto row_idx = 0u; row_idx < n_row; ++row_idx)
    for (auto col_idx = 0u; col_idx < n_col; ++col_idx)
      tmat[col_idx][row_idx] = mat[row_idx][col_idx];
}

std::istream &
operator>>(std::istream &is, Design &design) {
  std::string header_encoding;
  std::getline(is, header_encoding);

  std::istringstream header_is(header_encoding);
  std::string header_name;
  while (header_is >> header_name)
    design.factor_names.push_back(header_name);

  std::string row;
  while (std::getline(is, row)) {
    if (row.empty())
      continue;

    std::istringstream row_is(row);
    std::string token;
    row_is >> token;
    design.sample_names.push_back(token);

    std::vector<std::uint8_t> matrix_row;
    while (row_is >> token) {
      if (std::size(token) != 1 || (token != "0" && token != "1"))
        throw std::runtime_error("Must use binary factor levels:\n" + row);
      matrix_row.push_back(token == "1");
    }

    if (std::size(matrix_row) != design.n_factors())
      throw std::runtime_error(
        "each row must have as many columns as factors:\n" + row);

    design.matrix.push_back(matrix_row);
  }

  transpose(design.matrix, design.tmatrix);
  make_groups(design);
  assign_groups(design);

  return is;
}

[[nodiscard]] Design
Design::drop_factor(const std::uint32_t factor_idx) {
  // clang-format off
  Design design{
    factor_names,
    sample_names,
    matrix,
    tmatrix,
    groups,
    group_id,
  };
  // clang-format on
  design.factor_names.erase(std::begin(design.factor_names) + factor_idx);
  for (auto i = 0u; i < n_samples(); ++i)
    design.matrix[i].erase(std::begin(design.matrix[i]) + factor_idx);
  transpose(design.matrix, design.tmatrix);
  make_groups(design);
  assign_groups(design);
  return design;
}

void
Design::order_samples(const std::vector<std::string> &ordered_names) {
  // Build lookup: sample name -> original index
  const auto sample_to_index = [&] {
    std::unordered_map<std::string, std::uint32_t> sample_to_index;
    std::uint32_t index = 0;
    for (const auto &sample : sample_names)
      sample_to_index.emplace(sample, index++);
    return sample_to_index;
  }();

  std::vector<std::string> ord_sample_names;
  std::vector<std::vector<std::uint8_t>> ord_matrix;
  std::vector<std::uint32_t> ord_group_id;
  // factor names should not change
  // groups should not change
  // tmatrix will be changed after matrix

  const auto n_names = std::size(ordered_names);
  ord_sample_names.reserve(n_names);
  ord_matrix.reserve(n_names);
  ord_group_id.reserve(n_names);

  for (const auto &name : ordered_names) {
    const auto it = sample_to_index.find(name);
    if (it == std::cend(sample_to_index))
      throw std::runtime_error("Sample not found: " + name);

    const auto original_index = it->second;
    ord_sample_names.push_back(sample_names[original_index]);
    ord_matrix.push_back(matrix[original_index]);
    ord_group_id.push_back(group_id[original_index]);
  }
  sample_names = std::move(ord_sample_names);
  matrix = std::move(ord_matrix);
  group_id = std::move(ord_group_id);

  // update tmatrix using matrix
  transpose(matrix, tmatrix);
}

std::ostream &
operator<<(std::ostream &out, const Design &design) {
  static constexpr std::uint32_t max_samples_to_report = 20;
  const auto n_samples = design.n_samples();
  const auto n_factors = design.n_factors();
  if (n_samples <= max_samples_to_report) {
    for (std::size_t factor = 0; factor < n_factors; ++factor) {
      if (factor != 0)
        out << '\t';
      out << design.factor_names[factor];
    }
    out << '\n';

    for (std::size_t i = 0; i < n_samples; ++i) {
      out << design.sample_names[i];
      for (std::size_t j = 0; j < n_factors; ++j)
        out << '\t' << static_cast<std::uint32_t>(design.matrix[i][j]);
      out << '\n';
    }
  }

  // compute the number of samples per group
  const auto n_groups = design.n_groups();
  std::vector<std::uint32_t> n_samples_per_group(n_groups, 0);
  for (auto s_idx = 0u; s_idx < n_samples; ++s_idx)
    ++n_samples_per_group[design.group_id[s_idx]];

  out << "group_id | factor_levels (" << n_factors << " factors) | n_samples\n";
  for (std::size_t g_idx = 0; g_idx < n_groups; ++g_idx) {
    out << g_idx;
    for (std::size_t f_idx = 0; f_idx < n_factors; ++f_idx)
      out << '\t' << static_cast<std::uint32_t>(design.groups[g_idx][f_idx]);
    out << '\t' << n_samples_per_group[g_idx] << '\n';
  }

  return out;
}

void
ensure_sample_order(const std::string &table_filename, Design &design) {
  std::ifstream table_file(table_filename);
  if (!table_file)
    throw std::runtime_error("could not open file: " + table_filename);
  std::string header;
  std::getline(table_file, header);
  const auto sample_names = get_sample_names_from_header(header);
  design.order_samples(sample_names);
}

[[nodiscard]] std::vector<std::string>
get_sample_names_from_header(const std::string &header) {
  std::istringstream iss(header);
  std::string token;
  std::vector<std::string> sample_names;
  while (iss >> token)
    sample_names.push_back(token);
  return sample_names;
}

[[nodiscard]] bool
consistent_sample_names(const Design &design, const std::string &header) {
  std::istringstream iss(header);
  auto nm_itr(std::begin(design.sample_names));
  const auto nm_end(std::end(design.sample_names));
  std::string token;
  while (iss >> token && nm_itr != nm_end)
    if (token != *nm_itr++)
      return false;
  return true;
}

[[nodiscard]] Design
Design::read_design(const std::string &design_filename) {
  std::ifstream design_file(design_filename);
  if (!design_file)
    throw std::runtime_error("could not open file: " + design_filename);
  Design design;
  design_file >> design;
  return design;
}

[[nodiscard]] std::uint32_t
Design::get_test_factor_idx(const std::string &test_factor) const {
  const auto itr =
    std::find(std::cbegin(factor_names), std::cend(factor_names), test_factor);
  if (itr == std::cend(factor_names))
    throw std::runtime_error("factor not part of design: " + test_factor);
  return std::distance(std::cbegin(factor_names), itr);
}

[[nodiscard]] bool
Design::has_two_values(const std::size_t test_factor) const {
  const auto &tcol = tmatrix[test_factor];
  return std::any_of(std::cbegin(tcol), std::cend(tcol),
                     [&](const auto x) { return x != tcol[0]; });
}
