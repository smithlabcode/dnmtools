/* Copyright (C) 2018-2022 Andrew D. Smith
 *
 * Authors: Andrew D. Smith
 *
 * This is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * This software is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * for more details.
 */

#ifndef LEVELS_COUNTER_HPP
#define LEVELS_COUNTER_HPP

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
struct MSite;

struct LevelsCounter {

  // context is the dinucleotide or trinucleotide in which this
  // cytosine resides. The documentation for the counts and levels
  // tools provide more information about this context string.
  std::string context;

  // total_sites is the total number of sites summarized, all of which are
  // consistent with the context string. If the context string is empty or not
  // used, then any constraints on which sites are counted among this total is
  // defined separately.
  std::uint64_t total_sites{};

  // sites_covered is the number of reference genome sites analyzed
  // that have at least one mapped read covering them.
  std::uint64_t sites_covered{};

  // max_depth is the number of reads covering the reference genome site
  // having the most mapped reads covering it.
  std::uint64_t max_depth{};

  // mutations is the number of sites tagged as having been mutated in the
  // sample based on the mapped reads in comparisin with the underlying
  // reference genome site. This may be indicated in various ways and defined
  // in various ways.
  std::uint64_t mutations{};

  // total_c is the number of cytosines in mapped reads covering relevant
  // cytosines in the reference genome.
  std::uint64_t total_c{};

  // total_t is the number of thymines in mapped reads covering relevant
  // cytosines in the reference genome.
  std::uint64_t total_t{};

  // called_meth is the number of cytosines in the reference genome "called"
  // as being methylated using some statistical criteria.
  std::uint64_t called_meth{};

  // called_unmeth is the number of cytosines in the reference genome "called"
  // as being unmethylated using some statistical criteria.
  std::uint64_t called_unmeth{};

  // total_meth is the sum over all relevant cytosines of the methylation
  // level estimated at that cytosine using the mapped reads. This value
  // contributes to the unweighted mean methylation.
  double total_meth{};

  // coverage is equal to total_c plus total_t. These are the counts that
  // contribute towards estimates of the various methylation levels.
  [[nodiscard]] std::uint64_t
  coverage() const {
    return total_c + total_t;
  }

  // mean_depth is the average number of reads contributing to methylation
  // calls over all sites.
  [[nodiscard]] double
  mean_depth() const {
    return static_cast<double>(coverage()) / total_sites;
  }

  // mean_depth_covered is the average number of reads contributing to
  // methylation calls over all sites that are covered at least once.
  [[nodiscard]] double
  mean_depth_covered() const {
    return static_cast<double>(coverage()) / sites_covered;
  }

  // sites_covered_frac is the fraction of all sites that are covered at least
  // once.
  [[nodiscard]] double
  sites_covered_frac() const {
    return static_cast<double>(sites_covered) / total_sites;
  }

  // total_called is equal to called_meth plus called_unmeth
  [[nodiscard]] std::uint64_t
  total_called() const {
    return called_meth + called_unmeth;
  }

  // mean_meth_weighted is the unweighted mean methylation level. This
  // is the ratio of total_c divided by coverage. This value is always
  // between 0 and 1.
  [[nodiscard]] double
  mean_meth_weighted() const {
    return static_cast<double>(total_c) /
           std::max(coverage(), static_cast<std::uint64_t>(1));
  }

  // mean_meth is the unweighted mean methylation level. This is the ratio of
  // total_meth divided by sites_covered. Value is always between 0 and 1.
  [[nodiscard]] double
  mean_meth() const {
    return static_cast<double>(total_meth) /
           std::max(sites_covered, static_cast<std::uint64_t>(1));
  }

  // fractional_meth is the fraction of "called" sites that are called
  // methylated. It is the ratio of called_meth divided by total_called. This
  // value is always between 0 and 1.
  [[nodiscard]] double
  fractional_meth() const {
    return static_cast<double>(called_meth) /
           std::max(total_called(), static_cast<std::uint64_t>(1));
  }

  explicit LevelsCounter(const std::string &c) : context{c} {}

  LevelsCounter() = default;

  LevelsCounter &
  operator+=(const LevelsCounter &rhs);

  void
  update(const MSite &s);

  [[nodiscard]] std::string
  tostring() const;
  [[nodiscard]] std::string
  format_row() const;
  [[nodiscard]] static std::string
  format_header();

  static double alpha;
};

std::ostream &
operator<<(std::ostream &out, const LevelsCounter &cs);

std::istream &
operator>>(std::istream &in, LevelsCounter &cs);

#endif
