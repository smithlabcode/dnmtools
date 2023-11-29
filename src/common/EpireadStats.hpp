/* Copyright (C) 2011-2022 University of Southern California and
 *                    Andrew D. Smith and Fang Fang
 *
 * Authors: Fang Fang and Andrew D. Smith
 *
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 */

#ifndef EPIREAD_STATS
#define EPIREAD_STATS

#include <cstdint>
#include <vector>

#include "Epiread.hpp"

struct small_epiread {
  uint32_t pos{};
  std::string seq{};

  small_epiread(uint32_t p, std::string s): pos{p}, seq{s} {}

  uint32_t end() const { return pos + std::size(seq); }

  uint32_t length() const { return std::size(seq); }
};

double
log_likelihood(const small_epiread &r, const std::vector<double> &a);

void
fit_epiallele(const std::vector<small_epiread> &reads,
              const std::vector<double> &indicators, std::vector<double> &a);
double
fit_single_epiallele(const std::vector<small_epiread> &reads,
                     std::vector<double> &a);

double
log_likelihood(const small_epiread &r, const double z,
               const std::vector<double> &a1, const std::vector<double> &a2);
double
log_likelihood(const small_epiread &r, const std::vector<double> &a1,
               const std::vector<double> &a2);
double
log_likelihood(const std::vector<small_epiread> &reads,
               const std::vector<double> &indicators,
               const std::vector<double> &a1, const std::vector<double> &a2);

double
resolve_epialleles(const size_t max_itr,
                   const std::vector<small_epiread> &reads,
                   std::vector<double> &indicators, std::vector<double> &a1,
                   std::vector<double> &a2);

double
test_asm_lrt(const size_t max_itr, const bool crct_for_read_count,
             const double low_prob, const double high_prob,
             std::vector<small_epiread> &reads);

double
test_asm_bic(const size_t max_itr, const bool crct_for_read_count,
             const double low_prob, const double high_prob,
             std::vector<small_epiread> &reads);

struct EpireadStats {
  double test_asm(std::vector<small_epiread> &reads,
                  bool &is_significant) const {
    const double score = use_bic ? test_asm_bic(max_itr, crct_for_read_count,
                                                low_prob, high_prob, reads)
                                 : test_asm_lrt(max_itr, crct_for_read_count,
                                                low_prob, high_prob, reads);
    is_significant = use_bic ? score < 0.0 : score < critical_value;
    return score;
  }

  double low_prob{0.25};
  double high_prob{0.75};
  double critical_value{0.01};
  size_t max_itr{10};
  bool use_bic{false};
  bool crct_for_read_count{true};
};

#endif
