/* Copyright (C) 2025 Andrew D Smith
 *
 * Authors: Andrew D. Smith
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

#include <algorithm>
#include <charconv>
#include <cstdint>
#include <istream>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

double Regression::tolerance = 1e-4;
double Regression::stepsize = 0.01;
std::uint32_t Regression::max_iter = 250;

void
SiteProp::parse(const std::string &line) {
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
  mcounts mc1{};
  mc.clear();
  while (field_e != line_end) {
    // get the total count
    field_s = std::find_if(field_e + 1, line_end, not_sep);
    field_e = std::find_if(field_s + 1, line_end, is_sep);
    {
      const auto [ptr, ec] = std::from_chars(field_s, field_e, mc1.n_reads);
      failed = failed || (ptr != field_e);
    }

    field_s = std::find_if(field_e + 1, line_end, not_sep);
    field_e = std::find_if(field_s + 1, line_end, is_sep);
    {
      const auto [ptr, ec] = std::from_chars(field_s, field_e, mc1.n_meth);
      failed = failed || (ptr != field_e);
    }

    mc.push_back(mc1);
  }

  if (failed)
    throw std::runtime_error("failed to parse counts from:\n" + line);
}

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
Design::drop_factor(const std::size_t factor_idx) {
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

std::ostream &
operator<<(std::ostream &out, const Design &design) {
  static constexpr std::uint32_t max_samples_to_report = 200;
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
