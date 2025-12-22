/* Copyright (C) 2025 Andrew D Smith
 *
 * This is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * This is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this software; if not, write to the Free Software Foundation, Inc., 51
 * Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
 */

#ifndef INTERVAL_HPP_
#define INTERVAL_HPP_

#include <cstdint>
//  #include <format> // ADS: needs c++20
#include <fstream>
#include <iterator>  // std::size
#include <stdexcept>
#include <string>
#include <vector>

struct Interval {
  std::string chrom;
  std::uint32_t start{};
  std::uint32_t stop{};

  Interval() = default;
  Interval(const std::string &chrom, const std::uint32_t start,
           const std::uint32_t stop) : chrom{chrom}, start{start}, stop{stop} {}

  explicit Interval(const std::string &line) {
    if (!initialize(line.data(), line.data() + std::size(line)))
      throw std::runtime_error("bad interval line: " + line);
  }
  auto
  initialize(const char *, const char *) -> bool;

  [[nodiscard]] auto
  operator<(const Interval &rhs) const {
    return (chrom < rhs.chrom ||
            (chrom == rhs.chrom &&
             (start < rhs.start || (start == rhs.start && stop < rhs.stop))));
  }

  [[nodiscard]] auto
  operator==(const Interval &rhs) const {
    return chrom == rhs.chrom && start == rhs.start && stop < rhs.stop;
  }

  // auto
  // operator<=>(const Interval &) const = default;
};

inline auto
operator<<(std::ostream &os, const Interval &x) -> std::ostream & {
  return os << x.chrom << "\t" << x.start << "\t" << x.stop;
}

[[nodiscard]] inline auto
to_string(const Interval &x) -> std::string {
  return x.chrom + "\t" + std::to_string(x.start) + "\t" +
         std::to_string(x.stop);
}

// ADS: need to bump to c++20 for this
//
// template <> struct std::formatter<Interval> : std::formatter<std::string> {
//   auto
//   format(const Interval &i, format_context &ctx) const {
//     static constexpr auto fmt = "{}\t{}\t{}";
//     return std::formatter<std::string>::format(
//       std::format(fmt, i.chrom, i.start, i.stop), ctx);
//   }
// };

[[nodiscard]] inline auto
size(const Interval &x) {
  return x.stop > x.start ? x.stop - x.start : 0ul;
}

[[nodiscard]] auto
read_intervals(const std::string &intervals_file) -> std::vector<Interval>;

#endif  // INTERVAL_HPP_
