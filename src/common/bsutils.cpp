/* Copyright (C) 2018-2025 Andrew D. Smith
 *
 * Author: Andrew D. Smith
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

#include "bsutils.hpp"
#include "dnmtools_gaussinv.hpp"

#include <GenomicRegion.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

void
wilson_ci_for_binomial(const double alpha, const double n, const double p_hat,
                       double &lower, double &upper) {
  if (n <= 0.0) {  // protection
    lower = 0.0;
    upper = 1.0;
    return;
  }
  const double z = dnmt_gsl_cdf_ugaussian_Pinv(1 - alpha / 2);
  const double denom = 1 + z * z / n;
  const double first_term = p_hat + z * z / (2 * n);
  const double discriminant =
    std::max(0.0, p_hat * (1 - p_hat) / n + z * z / (4 * n * n));
  lower = std::max(0.0, (first_term - z * std::sqrt(discriminant)) / denom);
  upper = std::min(1.0, (first_term + z * std::sqrt(discriminant)) / denom);
}

void
adjust_region_ends(const std::vector<std::vector<GenomicRegion>> &clusters,
                   std::vector<GenomicRegion> &regions) {
  assert(std::size(clusters) == std::size(regions));
  for (std::size_t i = 0; i < std::size(regions); ++i) {
    std::size_t max_pos = regions[i].get_end();
    std::size_t min_pos = regions[i].get_start();
    for (std::size_t j = 0; j < std::size(clusters[i]); ++j) {
      max_pos = std::max(clusters[i][j].get_end(), max_pos);
      min_pos = std::min(clusters[i][j].get_start(), min_pos);
    }
    regions[i].set_end(max_pos);
    regions[i].set_start(min_pos);
  }
}

void
relative_sort(const std::vector<GenomicRegion> &mapped_locations,
              const std::vector<std::string> &names,
              std::vector<std::size_t> &lookup) {
  std::unordered_map<std::string, std::size_t> names_map;
  for (std::size_t i = 0; i < std::size(names); ++i)
    names_map[names[i]] = i;

  for (std::size_t i = 0; i < std::size(mapped_locations); ++i) {
    const auto j = names_map.find(mapped_locations[i].get_name());
    if (j == std::cend(names_map))
      throw std::runtime_error("read sequence not found for: " + names[i]);
    lookup.push_back(j->second);
  }
}
