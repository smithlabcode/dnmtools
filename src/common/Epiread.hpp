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

#ifndef EPIREAD
#define EPIREAD

#include "smithlab_utils.hpp"
#include <string>
#include <vector>

struct epiread {
  std::string chr{};
  std::size_t pos{};
  std::string seq{};
  epiread() = default;
  epiread(const std::string &line);
  epiread(const std::size_t p, const std::string &s) : pos(p), seq(s) {}
  epiread(const std::string &c, const std::size_t p, const std::string &s) :
    chr(c), pos(p), seq(s) {}
  epiread(const char *line, const std::size_t len);

  bool
  operator<(const epiread &other) const {
    return (chr < other.chr || (chr == other.chr && pos < other.pos));
  }
  std::size_t
  end() const {
    return pos + seq.length();
  }
  std::size_t
  length() const {
    return std::size(seq);
  }
  std::size_t
  size() const {
    return std::size(seq);
  }
};

std::istream &
operator>>(std::istream &in, epiread &er);
std::ostream &
operator<<(std::ostream &out, const epiread &er);

std::size_t
adjust_read_offsets(std::vector<epiread> &reads);

std::size_t
get_n_cpgs(const std::vector<epiread> &reads);

bool
validate_epiread_file(const std::string &filename);

#endif
