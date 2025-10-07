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

#include "LevelsCounter.hpp"
#include "bsutils.hpp"

#include <sstream>
#include <stdexcept>
#include <string>

using std::runtime_error;
using std::string;
using std::to_string;

void
LevelsCounter::update(const MSite &s) {
  if (s.is_mutated()) {
    ++mutations;
  }
  if (s.n_reads > 0) {
    ++sites_covered;
    max_depth = std::max(max_depth, static_cast<uint64_t>(s.n_reads));
    total_c += s.n_meth();
    total_t += s.n_reads - s.n_meth();
    total_meth += s.meth;
    double lower = 0.0, upper = 0.0;
    wilson_ci_for_binomial(alpha, s.n_reads, s.meth, lower, upper);
    called_meth += (lower > 0.5);
    called_unmeth += (upper < 0.5);
  }
  ++total_sites;
}

std::string
LevelsCounter::tostring() const {
  static const std::string indent = std::string(2, ' ');
  const bool good = (sites_covered != 0);
  std::ostringstream oss;
  // directly counted values
  oss << context + ":\n"
      << indent << "total_sites: " << total_sites << '\n'
      << indent << "sites_covered: " << sites_covered << '\n'
      << indent << "total_c: " << total_c << '\n'
      << indent << "total_t: " << total_t << '\n'
      << indent << "max_depth: " << max_depth << '\n'
      << indent << "mutations: " << mutations << '\n'
      << indent << "called_meth: " << called_meth << '\n'
      << indent << "called_unmeth: " << called_unmeth << '\n'
      << indent << "total_meth: " << total_meth << '\n';

  // derived values
  oss << indent << "coverage: " << coverage() << '\n'
      << indent << "sites_covered_fraction: " << sites_covered_fraction()
      << '\n'
      << indent << "mean_depth: " << mean_depth() << '\n'
      << indent << "mean_depth_covered: " << mean_depth_covered() << '\n'
      << indent << "mean_meth: " << (good ? std::to_string(mean_meth()) : "NA")
      << '\n'
      << indent << "mean_meth_weighted: "
      << (good ? std::to_string(mean_meth_weighted()) : "NA") << '\n'
      << indent << "fractional_meth: "
      << (good ? std::to_string(fractional_meth()) : "NA");
  return oss.str();
}

std::string
LevelsCounter::format_row() const {
  const bool good = (sites_covered > 0);
  std::ostringstream oss;
  // directly counted values
  // clang-format off
  oss << total_sites << '\t'
      << sites_covered << '\t'
      << total_c << '\t'
      << total_t << '\t'
      << max_depth << '\t'
      << mutations << '\t'
      << called_meth << '\t'
      << called_unmeth << '\t'
      << total_meth << '\t';
  // derived values
  oss << coverage() << '\t'
      << static_cast<double>(sites_covered)/total_sites << '\t'
      << static_cast<double>(coverage())/total_sites << '\t'
      << static_cast<double>(coverage())/sites_covered << '\t'
      << (good ? mean_meth() : 0.0)  << '\t'
      << (good ? mean_meth_weighted() : 0.0) << '\t'
      << (good ? fractional_meth() : 0.0);
  // clang-format on
  return oss.str();
}

std::string
LevelsCounter::format_header() {
  std::ostringstream oss;
  // directly counted values
  oss << "01. total_sites" << '\n'
      << "02. sites_covered" << '\n'
      << "03. total_c" << '\n'
      << "04. total_t" << '\n'
      << "05. max_depth" << '\n'
      << "06. mutations" << '\n'
      << "07. called_meth" << '\n'
      << "08. called_unmeth" << '\n'
      << "09. total_meth" << '\n';
  // derived values
  oss << "10. coverage" << '\n'
      << "11. sites_covered/total_sites" << '\n'
      << "12. coverage/total_sites" << '\n'
      << "13. coverage/sites_covered" << '\n'
      << "14. mean_meth" << '\n'
      << "15. mean_meth_weighted" << '\n'
      << "16. fractional_meth";
  return oss.str();
}

double LevelsCounter::alpha = 0.05;

std::ostream &
operator<<(std::ostream &out, const LevelsCounter &cs) {
  return out << cs.tostring();
}

static void
check_label(const std::string &observed, const std::string expected) {
  if (observed != expected)
    throw runtime_error("bad levels format [" + observed + "," + expected +
                        "]");
}

std::istream &
operator>>(std::istream &in, LevelsCounter &cs) {
  in >> cs.context;  // get the context
  cs.context = cs.context.substr(0, cs.context.find_first_of(":"));

  std::string label;
  in >> label >> cs.total_sites;  // the total sites
  check_label(label, "total_sites:");

  in >> label >> cs.sites_covered;  // the sites covered
  check_label(label, "sites_covered:");

  in >> label >> cs.total_c;  // the total c
  check_label(label, "total_c:");

  in >> label >> cs.total_t;  // the total t
  check_label(label, "total_t:");

  in >> label >> cs.max_depth;  // the max depth
  check_label(label, "max_depth:");

  in >> label >> cs.mutations;  // the number of mutations
  check_label(label, "mutations:");

  in >> label >> cs.called_meth;  // the number of sites called methylated
  check_label(label, "called_meth:");

  in >> label >> cs.called_unmeth;  // the number of sites called unmethylated
  check_label(label, "called_unmeth:");

  in >> label >> cs.total_meth;  // the mean aggregate
  check_label(label, "total_meth:");

  return in;
}

LevelsCounter &
LevelsCounter::operator+=(const LevelsCounter &rhs) {
  total_sites += rhs.total_sites;
  sites_covered += rhs.sites_covered;
  max_depth += rhs.max_depth;
  mutations += rhs.mutations;
  total_c += rhs.total_c;
  total_t += rhs.total_t;
  called_meth += rhs.called_meth;
  called_unmeth += rhs.called_unmeth;
  total_meth += rhs.total_meth;
  return *this;
}
