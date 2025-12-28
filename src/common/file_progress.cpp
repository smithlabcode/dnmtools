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

#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

// NOLINTBEGIN(*-narrowing-conversions)

[[nodiscard]] std::string
format_duration(const std::chrono::duration<double> elapsed) {
  static constexpr auto s_per_h = 3600;
  static constexpr auto s_per_m = 60;
  const double tot_s = elapsed.count();

  // break down into hours, minutes, seconds
  const std::uint32_t hours = tot_s / 3600;
  const std::uint32_t minutes = (static_cast<int>(tot_s) % s_per_h) / s_per_m;
  const double seconds = tot_s - (hours * s_per_h) - (minutes * s_per_m);

  std::ostringstream oss;
  // NOLINTBEGIN(*-avoid-magic-numbers)
  oss << std::setfill('0') << std::setw(2) << hours << ":" << std::setfill('0')
      << std::setw(2) << minutes << ":" << std::fixed << std::setprecision(2)
      << std::setw(5) << seconds;
  // NOLINTEND(*-avoid-magic-numbers)
  return oss.str();
}

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
