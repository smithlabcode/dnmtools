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

#ifndef INTERVAL6_HPP_
#define INTERVAL6_HPP_

#include <cstdint>
#include <format>
#include <iterator>  // std::size
#include <stdexcept>
#include <string>
#include <vector>

struct Interval6 {
  std::string chrom;
  std::uint32_t start{};
  std::uint32_t stop{};
  std::string name;
  double score{};
  char strand{};

  Interval6() = default;
  Interval6(const std::string &chrom, const std::uint32_t start,
            const std::uint32_t stop, const std::string &name,
            const double score, const char strand) :
    chrom{chrom},
    start{start}, stop{stop}, name{name}, score{score}, strand{strand} {}

  explicit Interval6(const std::string &line) {
    if (!initialize(line.data(), line.data() + std::size(line)))
      throw std::runtime_error("bad interval6 line: " + line);
  }
  auto
  initialize(const char *, const char *) -> bool;
  auto
  operator<=>(const Interval6 &) const = default;
};

template <> struct std::formatter<Interval6> : std::formatter<std::string> {
  auto
  format(const Interval6 &i, format_context &ctx) const {
    static constexpr auto fmt = "{}\t{}\t{}\t{}\t{:.6g}\t{}";
    return std::formatter<std::string>::format(
      std::format(fmt, i.chrom, i.start, i.stop, i.name, i.score, i.strand),
      ctx);
  }
};

[[nodiscard]] inline auto
size(const Interval6 &x) {
  return x.stop > x.start ? x.stop - x.start : 0ul;
}

[[nodiscard]] auto
read_intervals6(const std::string &intervals_file) -> std::vector<Interval6>;

#endif  // INTERVAL6_HPP_
