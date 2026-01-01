/* Copyright (C) 2025 Andrew D Smith
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

#include "file_progress.hpp"

#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

// NOLINTBEGIN(*-narrowing-conversions)

file_progress::file_progress(const std::string &filename) :
  one_thousand_over_filesize{
    static_cast<double>(one_thousand) /
    static_cast<double>(std::filesystem::file_size(filename))} {}

void
file_progress::operator()(
  std::ifstream &in) {  // cppcheck-suppress constParameterReference
  static constexpr auto ten = 10.0;
  const std::size_t curr_offset =
    in.eof() ? one_thousand : in.tellg() * one_thousand_over_filesize;
  if (curr_offset <= prev_offset)
    return;
  std::ios old_state(nullptr);
  old_state.copyfmt(std::cerr);
  std::cerr << "\r[progress: " << std::setw(5)  // NOLINT(*-avoid-magic-numbers)
            << std::fixed << std::setprecision(1) << (curr_offset / ten)
            << (curr_offset == one_thousand ? "%]\n" : "%]");
  std::cerr.copyfmt(old_state);
  prev_offset = curr_offset == one_thousand
                  ? std::numeric_limits<std::size_t>::max()
                  : curr_offset;
}

// NOLINTEND(*-narrowing-conversions)
