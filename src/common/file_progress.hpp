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

#ifndef SRC_COMMON_FILE_PROGRESS_HPP_
#define SRC_COMMON_FILE_PROGRESS_HPP_

#include <cstddef>
#include <fstream>
#include <string>

struct file_progress {
  static constexpr auto one_thousand{1000ul};
  double one_thousand_over_filesize{};
  std::size_t prev_offset{};
  explicit file_progress(const std::string &filename);
  void
  operator()(std::ifstream &in);  // cppcheck-suppress constParameterReference
};

#endif  // SRC_COMMON_FILE_PROGRESS_HPP_
