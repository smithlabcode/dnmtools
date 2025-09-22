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

#include "radmeth_summary_stats.hpp"

#include "radmeth_model.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

void
cache_dispersion_effects(Regression &reg, const std::vector<double> &phi) {
  // reg.distributed
  for (auto i = 0u; i < reg.n_groups(); ++i) {
    const auto max_k = reg.max_r_count[i];
    auto &cache = reg.disp_cache[i];
    const auto current_phi = phi[i];
    for (auto k = 0u; k < max_k; ++k)
      cache[k] = (k - 1.0) / (1.0 + current_phi * (k - 1.0));
  }
}

void
set_max_r_count(Regression &reg) {
  reg.disp_cache.resize(reg.n_groups());
  reg.max_r_count.resize(reg.n_groups());
  for (auto i = 0u; i < reg.n_groups(); ++i) {
    const auto &cumul = reg.cumul[i];
    reg.max_r_count[i] = std::size(cumul.r_counts);
    // ADS: avoid the realloc that can happen even for resize(smaller_size)
    if (reg.max_r_count[i] > std::size(reg.disp_cache[i]))
      reg.disp_cache[i].resize(reg.max_r_count[i]);
  }
}

void
cache_log1p_factors(Regression &reg, const std::vector<double> &phi) {
  for (auto i = 0u; i < std::size(phi); ++i) {
    const auto max_k = reg.max_r_count[i];
    auto &cache = reg.disp_cache[i];
    const auto current_phi = phi[i];
    for (auto k = 0u; k < max_k; ++k)
      cache[k] = std::log1p(current_phi * (k - 1.0));
  }
}

void
get_cumulative(const std::vector<std::uint32_t> &group_id,
               const std::uint32_t n_groups, const std::vector<mcounts> &mc,
               std::vector<cumul_counts> &cumul) {
  const auto n_cols = std::size(mc);
  cumul.clear();
  cumul.resize(n_groups);

  const auto comp_cumul = [&](auto get_value, auto get_vector) {
    // phase 1: determine max value for each group
    for (auto g_idx = 0u; g_idx < n_groups; ++g_idx) {
      std::uint32_t max_v{};
      for (auto c_idx = 0u; c_idx < n_cols; ++c_idx) {
        if (group_id[c_idx] == g_idx) {
          const auto val = get_value(mc[c_idx]);
          if (val > max_v)
            max_v = val;
        }
      }
      get_vector(cumul[g_idx]).resize(max_v, 0);
    }

    // phase 2: fill cumulative counts
    for (auto c_idx = 0u; c_idx < n_cols; ++c_idx) {
      const auto g_idx = group_id[c_idx];
      const auto val = get_value(mc[c_idx]);
      auto &vec = get_vector(cumul[g_idx]);
      for (auto i = 0u; i < val; ++i)
        ++vec[i];
    }
  };
  // call the lambda 3 times for m_counts, r_counts, d_counts
  comp_cumul(
    [](const mcounts &m) { return m.n_meth; },
    [](cumul_counts &c) -> std::vector<std::uint32_t> & { return c.m_counts; });

  comp_cumul(
    [](const mcounts &m) { return m.n_reads; },
    [](cumul_counts &c) -> std::vector<std::uint32_t> & { return c.r_counts; });

  comp_cumul(
    [](const mcounts &m) { return m.n_reads - m.n_meth; },
    [](cumul_counts &c) -> std::vector<std::uint32_t> & { return c.d_counts; });
}
