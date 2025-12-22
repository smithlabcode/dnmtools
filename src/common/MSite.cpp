/* Copyright (C) 2015-2025 Andrew D. Smith
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 2 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 */

#include "MSite.hpp"
#include "counts_header.hpp"

#include <bamxx.hpp>

#include <charconv>
#include <fstream>
#include <iterator>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>

bool MSite::no_extra_fields = true;

bool
MSite::initialize(const char *c, const char *c_end) {
  constexpr auto is_sep = [](const char x) { return x == ' ' || x == '\t'; };
  constexpr auto not_sep = [](const char x) { return x != ' ' && x != '\t'; };

  bool failed = false;

  // NOLINTBEGIN(*-pointer-arithmetic)
  auto field_s = c;
  auto field_e = std::find_if(field_s + 1, c_end, is_sep);
  if (field_e == c_end)
    failed = true;

  {
    const std::uint32_t d = std::distance(field_s, field_e);
    chrom = std::string{field_s, d};
  }

  field_s = std::find_if(field_e + 1, c_end, not_sep);
  field_e = std::find_if(field_s + 1, c_end, is_sep);
  failed = failed || (field_e == c_end);

  {
    const auto [ptr, ec] = std::from_chars(field_s, field_e, pos);
    failed = failed || (ptr == field_s);
  }

  field_s = std::find_if(field_e + 1, c_end, not_sep);
  field_e = std::find_if(field_s + 1, c_end, is_sep);
  failed = failed || (field_e != field_s + 1 || field_e == c_end);

  strand = *field_s;
  failed = failed || (strand != '-' && strand != '+');

  field_s = std::find_if(field_e + 1, c_end, not_sep);
  field_e = std::find_if(field_s + 1, c_end, is_sep);
  failed = failed || (field_e == c_end);

  {
    const std::uint32_t d = std::distance(field_s, field_e);
    context = std::string{field_s, d};
  }

  field_s = std::find_if(field_e + 1, c_end, not_sep);
  field_e = std::find_if(field_s + 1, c_end, is_sep);
  failed = failed || (field_e == c_end);

  {
#ifdef __APPLE__
    const int ret = std::sscanf(field_s, "%lf", &meth);
    failed = failed || (ret < 1);
#else
    const auto [ptr, ec] = std::from_chars(field_s, field_e, meth);
    failed = failed || (ptr == field_s);
#endif
  }

  field_s = std::find_if(field_e + 1, c_end, not_sep);

  {
    const auto [ptr, ec] = std::from_chars(field_s, c_end, n_reads);
    failed = failed || (no_extra_fields && ptr != c_end);
  }
  // ADS: the value below would work for a flag
  // pos = std::numeric_limits<decltype(pos)>::max();

  // NOLINTEND(*-pointer-arithmetic)

  return !failed;
}

MSite::MSite(const std::string &line) {
  // NOLINTNEXTLINE(*-pointer-arithmetic)
  if (!initialize(line.data(), line.data() + std::size(line)))
    throw std::runtime_error("bad count line: " + line);
}

MSite::MSite(const char *line, const int n) {
  // NOLINTNEXTLINE(*-pointer-arithmetic)
  if (!initialize(line, line + n))
    throw std::runtime_error("bad count line: " + std::string(line));
}

std::string
MSite::tostring() const {
  std::ostringstream oss;
  oss << chrom << '\t' << pos << '\t' << strand << '\t' << context << '\t'
      << meth << '\t' << n_reads;
  return oss.str();
}

size_t
distance(const MSite &a, const MSite &b) {
  return a.chrom == b.chrom ? std::max(a.pos, b.pos) - std::min(a.pos, b.pos)
                            : std::numeric_limits<size_t>::max();
}

static void
move_to_start_of_line(std::ifstream &in) {
  char next{};
  // move backwards by: one step forward, two steps back
  while (in.good() && in.get(next) && next != '\n') {
    in.unget();
    in.unget();
  }
  // ADS: need to reconsider this below
  if (in.bad())
    // hope this only happens when hitting the start of the file
    in.clear();
}

void
find_offset_for_msite(const std::string &chr, const size_t idx,
                      std::ifstream &site_in) {
  site_in.seekg(0, std::ios_base::beg);
  const auto begin_pos = site_in.tellg();
  site_in.seekg(0, std::ios_base::end);
  const auto end_pos = site_in.tellg();

  if (end_pos - begin_pos < 2)
    throw std::runtime_error("empty counts file");

  std::streamoff step_size = (end_pos - begin_pos) / 2;

  site_in.seekg(0, std::ios_base::beg);
  std::string low_chr;
  size_t low_idx{};
  if (!(site_in >> low_chr >> low_idx))
    throw std::runtime_error("failed search in counts file");

  // MAGIC: need the -2 here to get past the EOF and possibly a '\n'
  site_in.seekg(-2, std::ios_base::end);
  move_to_start_of_line(site_in);
  std::string high_chr;
  size_t high_idx{};
  if (!(site_in >> high_chr >> high_idx))
    throw std::runtime_error("failed search in counts file");

  std::streamoff pos = step_size;
  site_in.seekg(pos, std::ios_base::beg);
  move_to_start_of_line(site_in);
  std::streamoff prev_pos{};  // keep track of previous position in file

  // binary search inside sorted file on disk
  // iterate until step size is 0 or positions are identical
  while (step_size > 0 && prev_pos != pos) {
    // track (mid) position in file to make sure it keeps moving
    prev_pos = pos;

    std::string mid_chr;  // chromosome name at mid point
    size_t mid_idx{};     // position at mid point
    if (!(site_in >> mid_chr >> mid_idx))
      throw std::runtime_error("failed navigating inside file");

    // this check will never give a false indication of unsorted, but
    // might catch an unsorted file
    if (mid_chr < low_chr || mid_chr > high_chr) {
      // NOLINTBEGIN(performance-inefficient-string-concatenation)
      throw std::runtime_error("chromosomes unsorted inside file: "
                               "low=" +
                               low_chr + ",mid=" + mid_chr +
                               ",high=" + high_chr);
      // NOLINTEND(performance-inefficient-string-concatenation)
    }
    step_size /= 2;  // cut the range in half

    if (chr < mid_chr || (chr == mid_chr && idx <= mid_idx)) {
      // move to the left
      std::swap(high_chr, mid_chr);
      std::swap(high_idx, mid_idx);
      pos -= step_size;
    }
    else {
      // move to the left
      std::swap(low_chr, mid_chr);
      std::swap(low_idx, mid_idx);
      pos += step_size;
    }
    site_in.seekg(pos, std::ios_base::beg);
    move_to_start_of_line(site_in);
  }
}

void
find_offset_for_msite(
  const std::unordered_map<std::string, std::uint32_t> &chrom_order,
  const std::string &chr, const size_t idx, std::ifstream &site_in) {
  site_in.seekg(0, std::ios_base::beg);
  const auto begin_pos = site_in.tellg();
  site_in.seekg(0, std::ios_base::end);
  const auto end_pos = site_in.tellg();

  if (end_pos - begin_pos < 2)
    throw std::runtime_error("empty counts file");

  std::streamoff step_size = (end_pos - begin_pos) / 2;

  const auto chr_idx_itr = chrom_order.find(chr);
  if (chr_idx_itr == std::cend(chrom_order)) {
    site_in.seekg(0, std::ios_base::end);
    return;
  }
  const auto chr_idx = chr_idx_itr->second;

  site_in.seekg(0, std::ios_base::beg);
  std::string low_chr;
  size_t low_idx{};
  if (!(site_in >> low_chr >> low_idx))
    throw std::runtime_error("failed search in counts file");

  const auto low_chr_idx_itr = chrom_order.find(low_chr);
  if (low_chr_idx_itr == std::cend(chrom_order))
    throw std::runtime_error("inconsistent chromosome order info");
  auto low_chr_idx = low_chr_idx_itr->second;

  // MAGIC: need the -2 here to get past the EOF and hopefully a '\n'
  site_in.seekg(-2, std::ios_base::end);
  move_to_start_of_line(site_in);
  std::string high_chr;
  size_t high_idx{};
  if (!(site_in >> high_chr >> high_idx))
    throw std::runtime_error("failed search in counts file");

  const auto high_chr_idx_itr = chrom_order.find(high_chr);
  if (high_chr_idx_itr == std::cend(chrom_order))
    throw std::runtime_error("inconsistent chromosome order info");
  auto high_chr_idx = high_chr_idx_itr->second;

  std::streamoff pos = step_size;
  site_in.seekg(pos, std::ios_base::beg);
  move_to_start_of_line(site_in);
  std::streamoff prev_pos{};  // keep track of previous position in file

  // binary search inside sorted file on disk
  // iterate until step size is 0 or positions are identical
  while (step_size > 0 && prev_pos != pos) {
    // track (mid) position in file to make sure it keeps moving
    prev_pos = pos;

    std::string mid_chr;  // chromosome name at mid point
    size_t mid_idx{};     // position at mid point
    if (!(site_in >> mid_chr >> mid_idx))
      throw std::runtime_error("failed navigating inside file");

    const auto mid_chr_idx_itr = chrom_order.find(mid_chr);
    if (mid_chr_idx_itr == std::cend(chrom_order))
      throw std::runtime_error("inconsistent chromosome order info");
    const auto mid_chr_idx = mid_chr_idx_itr->second;

    // this check will never give a false indication of unsorted, but
    // might catch an unsorted file
    if (mid_chr_idx < low_chr_idx || mid_chr_idx > high_chr_idx) {
      // NOLINTBEGIN(performance-inefficient-string-concatenation)
      throw std::runtime_error("chromosomes unsorted inside file: "
                               "low=" +
                               std::to_string(low_chr_idx) +
                               ",mid=" + std::to_string(mid_chr_idx) +
                               ",high=" + std::to_string(high_chr_idx));
      // NOLINTEND(performance-inefficient-string-concatenation)
    }

    step_size /= 2;  // cut the range in half

    if (chr_idx < mid_chr_idx || (chr_idx == mid_chr_idx && idx <= mid_idx)) {
      // move to the left
      high_chr_idx = mid_chr_idx;
      // high_idx = mid_idx;
      pos -= step_size;
    }
    else {
      // move to the left
      low_chr_idx = mid_chr_idx;
      // low_idx = mid_idx;
      pos += step_size;
    }
    site_in.seekg(pos, std::ios_base::beg);
    move_to_start_of_line(site_in);
  }
}

bool
is_msite_line(const std::string &line) {
  std::istringstream iss(line);

  std::string chrom;
  if (!(iss >> chrom))
    return false;

  std::int64_t pos{};
  if (!(iss >> pos || pos < 0))
    return false;

  std::string strand;
  if (!(iss >> strand) || (strand.size() != 1) ||
      ((strand != "+") && (strand != "-")))
    return false;

  std::string context;
  // ADS: below, allowing any context
  // std::regex pattern("^C[pHWX][GH]$");
  if (!(iss >> context))
    return false;

  double level{};
  if (!(iss >> level) || level < 0.0 || level > 1.0)
    return false;

  std::int64_t n_reads{};
  if (!(iss >> n_reads || n_reads < 0))
    return false;

  std::string temp;
  if (MSite::no_extra_fields && iss >> temp)
    return false;
  else
    return true;
}

bool
is_msite_file(const std::string &file) {
  bamxx::bgzf_file in(file, "r");
  if (!in)
    throw std::runtime_error("cannot open file: " + file);

  std::string line;
  while (getline(in, line) && is_counts_header_line(line))
    ;

  return is_msite_line(line);
}
