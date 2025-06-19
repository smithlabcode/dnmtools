/*    Copyright (C) 2011-2022 University of Southern California and
 *                       Andrew D. Smith and Fang Fang
 *
 *    Authors: Fang Fang and Andrew D. Smith
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 */

#include <charconv>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>

#include "Epiread.hpp"

std::size_t
adjust_read_offsets(std::vector<epiread> &reads) {
  std::size_t first_read_offset = std::numeric_limits<std::size_t>::max();
  for (std::size_t i = 0; i < reads.size(); ++i)
    first_read_offset = std::min(reads[i].pos, first_read_offset);
  for (std::size_t i = 0; i < reads.size(); ++i)
    reads[i].pos -= first_read_offset;
  return first_read_offset;
}

std::size_t
get_n_cpgs(const std::vector<epiread> &reads) {
  std::size_t n_cpgs = 0;
  for (std::size_t i = 0; i < reads.size(); ++i)
    n_cpgs = std::max(n_cpgs, reads[i].end());
  return n_cpgs;
}

std::istream &
operator>>(std::istream &in, epiread &er) {
  std::string buffer;
  if (getline(in, buffer)) {
    std::istringstream is(buffer);
    if (!(is >> er.chr >> er.pos >> er.seq))
      throw std::runtime_error("malformed epiread line:\n" + buffer);
  }
  return in;
}

std::ostream &
operator<<(std::ostream &out, const epiread &er) {
  return out << er.chr << '\t' << er.pos << '\t' << er.seq;
}

bool
validate_epiread_file(const std::string &filename) {
  const std::size_t max_lines_to_validate = 10000;
  std::ifstream in(filename);
  if (!in)
    throw std::runtime_error("failed to open file: " + filename);

  std::string c, s, other;
  std::size_t p = 0;

  std::size_t n_lines = 0;
  std::string line;
  while (getline(in, line) && n_lines++ < max_lines_to_validate) {
    std::istringstream iss(line);
    if (!(iss >> c >> p >> s) || iss >> other)
      return false;
  }
  return true;
}

epiread::epiread(const std::string &line) {
  constexpr auto is_sep = [](const char x) { return x == ' ' || x == '\t'; };
  constexpr auto not_sep = [](const char x) { return x != ' ' && x != '\t'; };

  using std::distance;
  using std::find_if;
  using std::from_chars;

  bool failed = false;

  const auto c = line.data();
  const auto c_end = c + line.size();

  auto field_s = c;
  auto field_e = find_if(field_s + 1, c_end, is_sep);
  if (field_e == c_end)
    failed = true;

  chr = std::string{field_s, static_cast<uint32_t>(distance(field_s, field_e))};

  field_s = find_if(field_e + 1, c_end, not_sep);
  field_e = find_if(field_s + 1, c_end, is_sep);
  failed = failed || (field_e == c_end);

  const auto [ptr, ec] = from_chars(field_s, field_e, pos);
  failed = failed || (ptr == field_s);

  field_s = find_if(field_e + 1, c_end, not_sep);
  field_e = find_if(field_s + 1, c_end, is_sep);
  failed = failed || (field_e != c_end);

  seq = std::string{field_s, static_cast<uint32_t>(distance(field_s, field_e))};

  if (failed) {
    throw std::runtime_error("bad epiread line: " + line);
    // ADS: the value below would work for a flag
    // pos = std::numeric_limits<decltype(pos)>::max();
  }
}

epiread::epiread(const char *line, const std::size_t len) {
  constexpr auto is_sep = [](const char x) { return x == ' ' || x == '\t'; };
  constexpr auto not_sep = [](const char x) { return x != ' ' && x != '\t'; };

  bool failed = false;

  const auto c = line;
  const auto c_end = c + len;

  auto field_s = c;
  auto field_e = std::find_if(field_s + 1, c_end, is_sep);
  if (field_e == c_end)
    failed = true;

  chr = std::string{field_s,
                    static_cast<uint32_t>(std::distance(field_s, field_e))};

  field_s = std::find_if(field_e + 1, c_end, not_sep);
  field_e = std::find_if(field_s + 1, c_end, is_sep);
  failed = failed || (field_e == c_end);

  const auto [ptr, ec] = std::from_chars(field_s, field_e, pos);
  failed = failed || (ptr == field_s);

  field_s = std::find_if(field_e + 1, c_end, not_sep);
  field_e = std::find_if(field_s + 1, c_end, is_sep);
  failed = failed || (field_e != c_end);

  seq = std::string{field_s,
                    static_cast<uint32_t>(std::distance(field_s, field_e))};

  if (failed) {
    throw std::runtime_error("bad epiread line: " + std::string(line, len));
    // ADS: the value below would work for a flag
    // pos = std::numeric_limits<decltype(pos)>::max();
  }
}
