/* Copyright (C) 2025 Andrew D Smith
 *
 * Author: Andrew D. Smith
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

#ifndef RADMETH_UTILS_HPP
#define RADMETH_UTILS_HPP

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <string>

[[nodiscard]] std::string
format_duration(const std::chrono::duration<double> elapsed);

struct file_progress {
  double one_thousand_over_filesize{};
  std::size_t prev_offset{};
  explicit file_progress(const std::string &filename);
  void
  operator()(std::ifstream &in);  // cppcheck-suppress constParameterReference
};

[[nodiscard]] double
llr_test(const double null_loglik, const double full_loglik);

[[nodiscard]] inline double
overdispersion_factor(const std::uint32_t n_samples, const double dispersion) {
  return (n_samples - 1) / (dispersion + 1);
}

#endif  // RADMETH_UTILS_HPP
