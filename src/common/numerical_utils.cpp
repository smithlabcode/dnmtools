/*
 * Copyright (C) 2011-2022 University of Southern California
 *                    Andrew D Smith and Qiang Song
 * Author: Qiang Song and Andrew D. Smith
 *
 * This is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 */

#include "numerical_utils.hpp"

#include <algorithm>
#include <cmath>
#include <iterator>  // IWYU pragma: keep
#include <vector>

double
log_sum_log_vec(const std::vector<double> &vals, const size_t limit) {
  const auto x = std::max_element(
    std::cbegin(vals), std::cbegin(vals) + static_cast<std::ptrdiff_t>(limit));
  const double max_val = *x;
  const std::size_t max_idx = std::distance(std::cbegin(vals), x);
  double sum = 1.0;
  for (std::size_t i = 0; i < limit; ++i)
    if (i != max_idx)
      sum += std::exp(vals[i] - max_val);  // cppcheck-suppress useStlAlgorithm
  return max_val + std::log(sum);
}

double
log_sum_log(const std::vector<double>::const_iterator &begin,
            const std::vector<double>::const_iterator &end) {
  const auto max_itr = std::max_element(begin, end);
  const double max_val = *max_itr;
  double sum = 1.0;
  for (auto itr = begin; itr < end; ++itr)
    if (itr != max_itr)
      sum += std::exp(*itr - max_val);  // cppcheck-suppress useStlAlgorithm
  return max_val + std::log(sum);
}
