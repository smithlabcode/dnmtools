/* xcounts_utils: code for doing things with xcounts format and some
 * for counts format that is common to several tools.
 *
 * Copyright (C) 2023-2024 Andrew D. Smith
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

#ifndef XCOUNTS_UTILS_HPP
#define XCOUNTS_UTILS_HPP

#include <string>
#include <vector>
#include <cstdint>
#include <unordered_map>

struct xcounts_entry {
  uint64_t pos{};  // absolute position
  uint32_t n_meth{};
  uint32_t n_unmeth{};
};

std::unordered_map<std::string, std::vector<xcounts_entry>>
read_xcounts_by_chrom(const uint32_t n_threads, const std::string &xcounts_file);

bool
get_is_xcounts_file(const std::string &filename);

#endif
