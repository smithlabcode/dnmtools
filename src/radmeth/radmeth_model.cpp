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
double Regression::stepsize = 0.001;
std::uint32_t Regression::max_iter = 700;

void
SiteProportions::parse(const std::string &line) {
  const auto first_ws = line.find_first_of(" \t");

  // Parse the row name (must be like: "chr:position:strand:context")
  bool failed = false;

  const auto label_end = line.data() + first_ws;
  auto field_s = line.data();
  auto field_e = std::find(field_s + 1, label_end, ':');
  if (field_e == label_end)
    failed = true;

  {
    const std::uint32_t d = std::distance(field_s, field_e);
    chrom = std::string{field_s, d};
  }

  field_s = field_e + 1;
  field_e = std::find(field_s + 1, label_end, ':');
  failed = failed || (field_e == label_end);

  {
    const auto [ptr, ec] = std::from_chars(field_s, field_e, position);
    failed = failed || (ptr == field_s);
  }

  field_s = field_e + 1;
  field_e = std::find(field_s + 1, label_end, ':');
  failed = failed || (field_e != field_s + 1 || field_e == label_end);

  strand = *field_s;
  failed = failed || (strand != '-' && strand != '+');

  field_s = field_e + 1;
  field_e = std::find_if(field_s + 1, label_end,
                         [](const auto x) { return x == ' ' || x == '\t'; });
  failed = failed || (field_e != label_end);

  {
    const std::uint32_t d = std::distance(field_s, field_e);
    context = std::string{field_s, d};
  }

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

    design.matrix.push_back(std::vector<std::uint8_t>());
    swap(design.matrix.back(), matrix_row);
  }

  const auto n_row = std::size(design.matrix);
  const auto n_col = std::size(design.matrix.front());
  design.tmatrix.resize(n_col, std::vector<std::uint8_t>(n_row, 0.0));
  for (auto row_idx = 0u; row_idx < n_row; ++row_idx)
    for (auto col_idx = 0u; col_idx < n_col; ++col_idx)
      design.tmatrix[col_idx][row_idx] = design.matrix[row_idx][col_idx];

  return is;
}

std::ostream &
operator<<(std::ostream &out, const Design &design) {
  for (std::size_t factor = 0; factor < design.factor_names.size(); ++factor) {
    if (factor != 0)
      out << '\t';
    out << design.factor_names[factor];
  }
  out << '\n';

  for (std::size_t i = 0; i < design.n_samples(); ++i) {
    out << design.sample_names[i];
    for (std::size_t j = 0; j < design.n_factors(); ++j)
      out << "\t" << design.matrix[i][j];
    out << "\n";
  }
  return out;
}
