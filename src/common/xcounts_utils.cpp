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

#include "xcounts_utils.hpp"
#include "counts_header.hpp"
#include "dnmt_error.hpp"

#include <bamxx.hpp>

#include <htslib/sam.h>

#include <cctype>
#include <charconv>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <system_error>
#include <unordered_map>
#include <vector>

using std::runtime_error;
using std::string;
using std::to_string;
using std::unordered_map;
using std::vector;

using bamxx::bgzf_file;

// careful: this could get big
unordered_map<string, vector<xcounts_entry>>
read_xcounts_by_chrom(const std::int32_t n_threads,
                      const string &xcounts_file) {
  bamxx::bam_tpool tp(n_threads);
  bamxx::bgzf_file in(xcounts_file, "r");
  if (!in)
    throw runtime_error("failed to open input file");

  // set the threads for the input file decompression
  if (n_threads > 1 && in.is_bgzf())
    tp.set_io(in);

  kstring_t line{0, 0, nullptr};
  string chrom_name;
  uint64_t pos = 0;

  unordered_map<string, vector<xcounts_entry>> sites_by_chrom;

  vector<xcounts_entry> curr_chrom;

  while (bamxx::getline(in, line)) {
    if (is_counts_header_line(line.s))
      continue;  // ADS: early loop exit

    // check if we have a chrom line
    if (!std::isdigit(line.s[0])) {  // NOLINT(*-pointer-arithmetic)
      if (!chrom_name.empty()) {
        sites_by_chrom.insert({chrom_name, curr_chrom});
        curr_chrom.clear();
      }
      chrom_name = string{line.s};
      pos = 0;
      continue;
    }

    std::uint32_t pos_step{};
    std::uint32_t n_meth{};
    std::uint32_t n_unmeth{};
    // NOLINTBEGIN(*-pointer-arithmetic)
    const auto end_line = line.s + line.l;
    auto res = std::from_chars(line.s, end_line, pos_step);
    res = std::from_chars(res.ptr + 1, end_line, n_meth);
    res = std::from_chars(res.ptr + 1, end_line, n_unmeth);
    // NOLINTEND(*-pointer-arithmetic)

    const auto curr_pos = pos + pos_step;

    curr_chrom.push_back({curr_pos, n_meth, n_unmeth});

    pos = curr_pos;
  }

  if (!chrom_name.empty())
    sites_by_chrom.insert({chrom_name, curr_chrom});

  ks_free(&line);
  return sites_by_chrom;
}

bool
get_is_xcounts_file(const std::string &filename) {
  static constexpr auto max_lines_to_check = 1000ul;
  bamxx::bgzf_file in(filename, "r");
  if (!in)
    throw dnmt_error{"failed to open input file: " + filename};

  kstring_t line{0, 0, nullptr};

  bool found_header = false;
  bool found_chrom = false;

  auto line_count = 0ul;
  while (line_count++ < max_lines_to_check && bamxx::getline(in, line)) {
    if (is_counts_header_line(line.s)) {
      found_header = true;
    }
    // check if we have a chrom line
    else if (!std::isdigit(line.s[0])) {  // NOLINT(*-pointer-arithmetic)
      if (!found_header) {
        ks_free(&line);
        return false;
      }
      found_chrom = true;
    }
    else {
      if (!found_chrom) {
        ks_free(&line);
        return false;
      }
      std::int64_t pos_step = 0, n_meth = 0, n_unmeth = 0;
      // NOLINTBEGIN(*-pointer-arithmetic)
      const auto end_line = line.s + line.l;
      auto res = std::from_chars(line.s, end_line, pos_step);
      if (res.ec != std::errc()) {
        ks_free(&line);
        return false;
      }
      res = std::from_chars(res.ptr + 1, end_line, n_meth);
      if (res.ec != std::errc()) {
        ks_free(&line);
        return false;
      }
      res = std::from_chars(res.ptr + 1, end_line, n_unmeth);
      if (res.ec != std::errc()) {
        ks_free(&line);
        return false;
      }
      // NOLINTEND(*-pointer-arithmetic)
    }
  }
  ks_free(&line);
  return true;
}
