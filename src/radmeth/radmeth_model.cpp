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
#include <string>
#include <vector>

double Regression::tolerance = 1e-4;
double Regression::stepsize = 0.01;
std::uint32_t Regression::max_iter = 250;

void
Regression::parse(const std::string &line) {
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
