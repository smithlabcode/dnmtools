/* Copyright (C) 2018-2022 Andrew D. Smith
 *
 * Authors: Andrew D. Smith
 *
 * This is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 */

#ifndef LEVELS_COUNTER_HPP
#define LEVELS_COUNTER_HPP

#include "MSite.hpp"

#include <iostream>
#include <string>

struct LevelsCounter {
  std::string context;
  uint64_t total_sites{};
  uint64_t sites_covered{};
  uint64_t max_depth{};
  uint64_t mutations{};
  uint64_t total_c{}, total_t{};
  uint64_t called_meth{}, called_unmeth{};
  double total_meth{};
  LevelsCounter(const std::string &c) : context{c} {}

  LevelsCounter() = default;

  LevelsCounter &operator+=(const LevelsCounter &rhs);

  void update(const MSite &s);

  uint64_t coverage() const {return total_c + total_t;}
  uint64_t total_called() const {return called_meth + called_unmeth;}

  double mean_meth_weighted() const {
    return static_cast<double>(total_c)/std::max(coverage(), 1ul);
  }
  double fractional_meth() const {
    return static_cast<double>(called_meth)/std::max(total_called(), 1ul);
  }
  double mean_meth() const {
    return static_cast<double>(total_meth)/std::max(sites_covered, 1ul);
  }

  std::string tostring() const;
  std::string format_row() const;
  static std::string format_header();

  static double alpha;
};

std::ostream &
operator<<(std::ostream &out, const LevelsCounter &cs);

std::istream &
operator>>(std::istream &in, LevelsCounter &cs);

#endif
