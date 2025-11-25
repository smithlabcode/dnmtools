/*    Copyright (C) 2014-2022 University of Southern California and
 *                       Andrew D. Smith and Fang Fang and Benjamin Decato
 *
 *    Authors: Fang Fang and Benjamin Decato and Andrew D. Smith
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

#include "EpireadStats.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <limits>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

#include <gsl/gsl_cdf.h>

using std::isfinite;
using std::max;
using std::min;
using std::string;
using std::vector;

template <typename T> using num_lim = std::numeric_limits<T>;
using epi_r = small_epiread;

static const double PSEUDOCOUNT = 1e-10;

static inline uint32_t
adjust_read_offsets(vector<epi_r> &reads) {
  const auto first_read_offset = std::accumulate(
    std::cbegin(reads), std::cend(reads), num_lim<std::uint32_t>::max(),
    [](const std::uint32_t a, const auto &r) { return std::min(a, r.pos); });
  for (auto &r : reads)
    r.pos -= first_read_offset;
  return first_read_offset;
}

static inline uint32_t
get_n_cpgs(const vector<epi_r> &reads) {
  return std::accumulate(
    std::cbegin(reads), std::cend(reads), 0u,
    [](const std::uint32_t a, const auto &r) { return std::max(a, r.end()); });
}

inline double
log_likelihood(const epi_r &r, const vector<double> &a) {
  double ll = 0.0;
  for (size_t i = 0; i < r.seq.length(); ++i)
    if (r.seq[i] == 'C' || r.seq[i] == 'T')
      ll += log(r.seq[i] == 'C' ? a[r.pos + i] : (1.0 - a[r.pos + i]));
  assert(std::isfinite(ll));
  return ll;
}

// inline double
// log_likelihood(const epi_r &r, const vector<double> &a) {
//   double ll = 0.0;
//   auto a_itr = cbegin(a) + r.pos;
//   for (auto s : r.seq) {
//     if (s != 'N')
//       ll += log((s == 'C') ? *a_itr : 1.0 - *a_itr);
//     ++a_itr;
//   }
//   return ll;
// }

static inline std::pair<double, double>
log_likelihood_pair(const epi_r &r, const vector<double> &a1,
                    const vector<double> &a2) {
  double ll1 = 0.0, ll2 = 0.0;
  // auto a1_itr = cbegin(a1) + r.pos;
  // auto a2_itr = cbegin(a2) + r.pos;
  // for (auto s : r.seq) {
  for (size_t i = 0u; i < r.seq.length(); ++i) {
    if (r.seq[i] == 'C') {
      ll1 += log(a1[r.pos + i]);
      ll2 += log(a2[r.pos + i]);
    }
    if (r.seq[i] == 'T') {
      ll1 += log(1.0 - a1[r.pos + i]);
      ll2 += log(1.0 - a2[r.pos + i]);
    }
    // ++a1_itr;
    // ++a2_itr;
  }
  return {ll1, ll2};
}

inline double
log_likelihood(const epi_r &r, const double log_mixing1,
               const double log_mixing2, const vector<double> &a1,
               const vector<double> &a2) {
  auto [ll1, ll2] = log_likelihood_pair(r, a1, a2);
  return log(exp(log_mixing1 + ll1) + exp(log_mixing2 + ll2));
}

inline double
log_likelihood(const vector<epi_r> &reads, const double log_mixing1,
               const double log_mixing2, const vector<double> &a1,
               const vector<double> &a2) {
  return std::transform_reduce(
    cbegin(reads), cend(reads), 0.0, std::plus<>(), [&](const epi_r &r) {
      return log_likelihood(r, log_mixing1, log_mixing2, a1, a2);
    });
}

static double
expectation_step(const vector<epi_r> &reads, const double mixing,
                 const vector<double> &a1, const vector<double> &a2,
                 vector<double> &indicators) {
  const double log_mixing1 = log(mixing);
  const double log_mixing2 = log(1.0 - mixing);
  assert(isfinite(log_mixing1) && isfinite(log_mixing2));

  double score = 0.0;
  auto ind_itr = begin(indicators);
  for (const auto &r : reads) {
    // for (uint32_t i = 0; i < reads.size(); ++i) {
    const double ll1 = log_mixing1 + log_likelihood(r, a1);
    const double ll2 = log_mixing2 + log_likelihood(r, a2);
    // assert(isfinite(ll1) && isfinite(ll2));
    const double log_denom = log(exp(ll1) + exp(ll2));
    score += log_denom;
    *ind_itr++ = exp(ll1 - log_denom);
    // assert(isfinite(log_denom) && isfinite(indicators[i]));
  }
  return score;
}

template <const bool inverse>
void
fit_epiallele(const double pseudo, const vector<epi_r> &reads,
              vector<double>::const_iterator indic_itr, vector<double> &a) {
  vector<double> total(size(a), 2 * pseudo);
  auto t_beg = begin(total);
  auto a_beg = begin(a);
  std::fill_n(a_beg, std::size(a), pseudo);
  for (const auto &r : reads) {
    const double weight = inverse ? 1.0 - *indic_itr++ : *indic_itr++;
    auto m_itr = a_beg + r.pos;
    auto t_itr = t_beg + r.pos;
    for (auto s : r.seq) {
      *m_itr++ += weight * (s == 'C');
      *t_itr++ += weight * (s != 'N');
    }
  }
  while (t_beg != cend(total))
    *a_beg++ /= *t_beg++;
}

void
fit_epiallele(const double pseudo, const vector<epi_r> &reads,
              vector<double> &a) {
  const uint32_t n_cpgs = a.size();
  vector<double> total(n_cpgs, 2 * pseudo);
  auto t_beg = begin(total);
  auto a_beg = begin(a);
  std::fill_n(a_beg, n_cpgs, pseudo);
  for (const auto &r : reads) {
    const uint32_t start = r.pos;
    auto m_itr = a_beg + start;
    auto t_itr = t_beg + start;
    for (const auto s : r.seq) {
      *m_itr++ += (s == 'C');
      *t_itr++ += (s != 'N');
    }
  }
  while (t_beg != cend(total))
    *a_beg++ /= *t_beg++;
}

static inline void
rescale_indicators(
  const double mixing,
  vector<double> &indic) {  // cppcheck-suppress constParameterReference
  const double ratio = std::reduce(std::cbegin(indic), std::cend(indic), 0.0) /
                       static_cast<double>(std::size(indic));
  if (mixing < ratio) {
    const double adjustment = mixing / ratio;
    std::transform(std::cbegin(indic), std::cend(indic), std::begin(indic),
                   [&](const auto x) { return x * adjustment; });
  }
  else {
    const double adjustment = mixing / (1.0 - ratio);
    std::transform(std::cbegin(indic), std::cend(indic), std::begin(indic),
                   [&](const auto x) { return 1.0 - (1.0 - x) * adjustment; });
  }
}

static double
expectation_maximization(const size_t max_itr, const vector<epi_r> &reads,
                         const double &mixing, vector<double> &indicators,
                         vector<double> &a1, vector<double> &a2) {
  static constexpr double EPIREAD_STATS_TOLERANCE = 1e-10;

  double prev_score = -num_lim<double>::max();
  for (auto i = 0u; i < max_itr; ++i) {
    const double score = expectation_step(reads, mixing, a1, a2, indicators);
    rescale_indicators(mixing, indicators);

    if ((prev_score - score) / prev_score < EPIREAD_STATS_TOLERANCE)
      break;

    // maximization_step(reads, indicators, a1, a2);

    // NOLINTBEGIN(*-avoid-magic-numbers)
    fit_epiallele<false>(0.5 * PSEUDOCOUNT, reads, cbegin(indicators), a1);
    fit_epiallele<true>(0.5 * PSEUDOCOUNT, reads, cbegin(indicators), a2);
    // NOLINTEND(*-avoid-magic-numbers)

    prev_score = score;
  }
  return prev_score;
}

double
resolve_epialleles(const size_t max_itr, const vector<epi_r> &reads,
                   const double &mixing, vector<double> &indicators,
                   vector<double> &a1, vector<double> &a2) {
  indicators.clear();
  indicators.resize(reads.size(), 0.0);
  for (size_t i = 0; i < reads.size(); ++i) {
    const double l1 = log_likelihood(reads[i], a1);
    const double l2 = log_likelihood(reads[i], a2);
    indicators[i] = exp(l1 - log(exp(l1) + exp(l2)));
  }

  return expectation_maximization(max_itr, reads, mixing, indicators, a1, a2);
}

double
fit_single_epiallele(const vector<epi_r> &reads, vector<double> &a) {
  // assert(reads.size() > 0);
  fit_epiallele(PSEUDOCOUNT, reads, a);

  double score = 0.0;
  for (const auto &r : reads) {
    score += log_likelihood(r, a);  // cppcheck-suppress useStlAlgorithm
    // assert(isfinite(score));
  }
  return score;
}

void
compute_model_likelihoods(double &single_score, double &pair_score,
                          const size_t &max_itr, const double &low_prob,
                          const double &high_prob, const size_t n_cpgs,
                          const vector<epi_r> &reads) {
  static constexpr double mixing = 0.5;

  // try a single epi-allele and compute its log likelihood
  vector<double> a0(n_cpgs, 0.5);  // NOLINT(*-magic-numbers)
  single_score = fit_single_epiallele(reads, a0);

  // initialize the pair epi-alleles and indicators, and do the actual
  // computation to infer alleles, compute its log likelihood
  vector<double> a1(n_cpgs, low_prob), a2(n_cpgs, high_prob), indicators;
  resolve_epialleles(max_itr, reads, mixing, indicators, a1, a2);

  const double log_mixing1 = log(mixing);
  const double log_mixing2 = log(1.0 - mixing);

  pair_score = transform_reduce(
    cbegin(reads), cend(reads), 0.0, std::plus<>(), [&](const epi_r &r) {
      return log_likelihood(r, log_mixing1, log_mixing2, a1, a2);
    });
}

double
test_asm_lrt(const size_t max_itr, const bool correct_for_read_count,
             const double low_prob, const double high_prob,
             vector<epi_r> &reads) {
  // NOLINTBEGIN(*-avoid-magic-numbers)
  double single_score = num_lim<double>::min();
  double pair_score = num_lim<double>::min();
  const auto first_read_offset = adjust_read_offsets(reads);
  const auto n_cpgs = get_n_cpgs(reads);

  compute_model_likelihoods(single_score, pair_score, max_itr, low_prob,
                            high_prob, n_cpgs, reads);

  for (auto &read : reads)
    read.pos += first_read_offset;

  // degrees of freedom = 2*n_cpgs for two-allele model
  // minus n_cpgs for one-allele model
  const size_t df = n_cpgs;

  // correction for numbers of reads
  if (correct_for_read_count)
    pair_score += static_cast<double>(std::size(reads)) * std::log(0.5);

  const double llr_stat = -2.0 * (single_score - pair_score);
  const double p_value =
    1.0 - gsl_cdf_chisq_P(llr_stat, static_cast<double>(df));
  // NOLINTEND(*-avoid-magic-numbers)

  return p_value;
}

double
test_asm_bic(const size_t max_itr, const bool correct_for_read_count,
             const double low_prob, const double high_prob,
             vector<epi_r> &reads) {
  // NOLINTBEGIN(*-avoid-magic-numbers)
  double single_score = num_lim<double>::min();
  double pair_score = num_lim<double>::min();
  const auto first_read_offset = adjust_read_offsets(reads);
  const auto n_cpgs = get_n_cpgs(reads);

  compute_model_likelihoods(single_score, pair_score, max_itr, low_prob,
                            high_prob, n_cpgs, reads);

  // correction for numbers of reads
  if (correct_for_read_count)
    pair_score += static_cast<double>(std::size(reads)) * log(0.5);

  for (auto &read : reads)
    read.pos += first_read_offset;

  const double n_reads = static_cast<double>(std::size(reads));
  const auto n_cpgs_log_n_reads = n_cpgs * std::log(n_reads);

  // compute bic scores and compare
  const double bic_single = n_cpgs_log_n_reads - 2 * single_score;
  const double bic_pair = 2 * n_cpgs_log_n_reads - 2 * pair_score;
  return bic_pair - bic_single;
  // NOLINTEND(*-avoid-magic-numbers)
}
