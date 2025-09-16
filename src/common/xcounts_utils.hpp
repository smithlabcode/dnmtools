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

#include <cstdint>
#include <ostream>
#include <string>
#include <unordered_map>
#include <vector>

struct xcounts_entry {
  std::uint64_t pos{};  // absolute position
  std::uint32_t n_meth{};
  std::uint32_t n_unmeth{};

  [[nodiscard]] std::uint32_t
  n_reads() const {
    return n_meth + n_unmeth;
  }

  [[nodiscard]] double
  frac() const {
    return static_cast<double>(n_meth) / n_reads();
  }
};

inline std::ostream &
operator<<(std::ostream &o, const xcounts_entry &e) {
  return o << e.pos << '\t' << e.n_meth << '\t' << e.n_unmeth;
}

std::unordered_map<std::string, std::vector<xcounts_entry>>
read_xcounts_by_chrom(const std::uint32_t n_threads,
                      const std::string &xcounts_file);

bool
get_is_xcounts_file(const std::string &filename);

#endif
