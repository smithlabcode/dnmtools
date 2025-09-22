/* Copyright (C) 2025 Andrew D Smith
 *
 * Authors: Andrew D. Smith
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 */

#ifndef RADMETH_SUMMARY_STATS_HPP
#define RADMETH_SUMMARY_STATS_HPP

#include "radmeth_model.hpp"

#include <cstdint>
#include <istream>
#include <ostream>
#include <string>
#include <vector>

void
cache_dispersion_effects(Regression &reg, const std::vector<double> &phis);

void
cache_log1p_factors(Regression &reg, const std::vector<double> &phis);

void
set_max_r_count(Regression &reg);

void
get_cumulative(const std::vector<std::uint32_t> &group_id,
               const std::uint32_t n_groups, const std::vector<mcounts> &mc,
               std::vector<cumul_counts> &cumul);

#endif  // RADMETH_SUMMARY_STATS_HPP
