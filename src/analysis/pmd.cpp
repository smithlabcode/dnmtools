/* Copyright (C) 2018-2025 Andrew D Smith and Benjamin Decato
 *
 * This is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * This is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 */

[[maybe_unused]] static constexpr auto description = R"(
Identify PMDs in methylomes. Methylation must be provided in the methcounts
file format (chrom, position, strand, context, methylation, reads). See the
methcounts documentation for details. This program assumes only data at CpG
sites and that strands are collapsed so only the positive site appears in the
file, but reads counts are from both strands.
)";

#include "Interval.hpp"
#include "Interval6.hpp"
#include "MSite.hpp"
#include "TwoStateHMM_PMD.hpp"
#include "bsutils.hpp"
#include "counts_header.hpp"

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"

#include <bamxx.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <new>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

// NOLINTBEGIN(*-narrowing-conversions)

struct pmd_summary {
  explicit pmd_summary(const std::vector<Interval6> &pmds) :
    pmd_count{std::size(pmds)},
    pmd_total_size{std::accumulate(
      std::cbegin(pmds), std::cend(pmds), 0ul,
      [](const std::uint64_t t, const Interval6 &p) { return t + size(p); })} {
    pmd_mean_size = static_cast<double>(pmd_total_size) /
                    std::max(pmd_count, static_cast<std::uint64_t>(1));
  }
  // pmd_count is the number of identified PMDs.
  std::uint64_t pmd_count{};
  // total_pmd_size is the sum of the sizes of the identified PMDs
  std::uint64_t pmd_total_size{};
  // mean_pmd_size is the mean size of the identified PMDs
  double pmd_mean_size{};
};

[[nodiscard]] static inline auto
to_string(const pmd_summary &s) -> std::string {
  std::ostringstream oss;
  oss << "pmd_count: " << s.pmd_count << '\n'
      << "pmd_total_size: " << s.pmd_total_size << '\n'
      << "pmd_mean_size: " << std::fixed << s.pmd_mean_size;
  return oss.str();
}

static void
get_adjacent_distances(const std::vector<Interval6> &pmds,
                       std::vector<std::size_t> &dists) {
  for (std::size_t i = 1; i < std::size(pmds); ++i)
    if (pmds[i].chrom == pmds[i - 1].chrom)
      dists.push_back(pmds[i].start - pmds[i - 1].stop);
}

static bool
precedes(const std::string &chrom, const std::size_t position,
         const Interval6 &r) {
  return chrom < r.chrom || (chrom == r.chrom && position < r.start);
}

static bool
succeeds(const std::string &chrom, const std::size_t position,
         const Interval6 &r) {
  return r.chrom < chrom || (chrom == r.chrom && r.stop <= position);
}

static void
merge_nearby_pmd(const std::size_t max_merge_dist,
                 std::vector<Interval6> &pmds) {
  std::size_t j = 0;
  for (std::size_t i = 1; i < std::size(pmds); ++i) {
    if (pmds[j].chrom == pmds[i].chrom &&
        pmds[i].start - pmds[j].stop <= max_merge_dist) {
      pmds[j].stop = pmds[i].stop;
      const std::string combined_name(pmds[j].name + pmds[i].name);
      pmds[j].name = combined_name;
    }
    else
      pmds[++j] = pmds[i];
  }
  pmds.resize(j);
}

[[nodiscard]] static inline double
lnbeta(const double a, const double b) {
  return std::lgamma(a) + std::lgamma(b) - std::lgamma(a + b);
}

[[nodiscard]] static inline double
beta_log_likelihood(const double alpha, const double beta, const double p) {
  return (alpha - 1.0) * std::log(p) + (beta - 1.0) * std::log(1.0 - p) -
         lnbeta(alpha, beta);
}

[[nodiscard]] static inline double
beta_max_likelihood(const double fg_alpha, const double fg_beta,
                    const double bg_alpha, const double bg_beta,
                    const double p_low, const double p_hi) {
  return beta_log_likelihood(fg_alpha, fg_beta, p_low) +
         beta_log_likelihood(bg_alpha, bg_beta, p_hi);
}

static std::size_t
find_best_bound(const bool IS_RIGHT_BOUNDARY,
                const std::map<std::size_t, std::pair<std::size_t, std::size_t>>
                  &pos_meth_tot,
                const std::vector<double> &fg_alpha,
                const std::vector<double> &fg_beta,
                const std::vector<double> &bg_alpha,
                const std::vector<double> &bg_beta) {
  std::vector<std::pair<std::size_t, std::size_t>> meth_tot;
  std::vector<std::size_t> positions;
  for (const auto &it : pos_meth_tot) {
    positions.push_back(it.first);
    meth_tot.push_back(it.second);
  }

  std::vector<std::size_t> cumu_left_meth(std::size(meth_tot), 0);
  std::vector<std::size_t> cumu_left_tot(std::size(meth_tot), 0);
  std::vector<std::size_t> cumu_right_meth(std::size(meth_tot), 0);
  std::vector<std::size_t> cumu_right_tot(std::size(meth_tot), 0);
  if (std::size(meth_tot) > 0)
    for (std::size_t i = 1; i + 1 < std::size(meth_tot); ++i) {
      const std::size_t j = std::size(meth_tot) - 1 - i;
      cumu_left_meth[i] = cumu_left_meth[i - 1] + meth_tot[i - 1].first;
      cumu_left_tot[i] = cumu_left_tot[i - 1] + meth_tot[i - 1].second;
      cumu_right_meth[j] = cumu_right_meth[j + 1] + meth_tot[j + 1].first;
      cumu_right_tot[j] = cumu_right_tot[j + 1] + meth_tot[j + 1].second;
    }

  std::size_t best_idx = 0;
  double best_score = -std::numeric_limits<double>::max();
  if (std::size(meth_tot) > 0)
    for (std::size_t i = 1; i + 1 < std::size(meth_tot); ++i) {
      std::size_t N_low{};
      std::size_t k_low{};
      std::size_t N_hi{};
      std::size_t k_hi{};
      if (!IS_RIGHT_BOUNDARY) {
        N_low = cumu_right_tot[i] + meth_tot[i].second;
        k_low = cumu_right_meth[i] + meth_tot[i].first;
        N_hi = cumu_left_tot[i];
        k_hi = cumu_left_meth[i];
      }
      else {
        N_low = cumu_left_tot[i] + meth_tot[i].second;
        k_low = cumu_left_meth[i] + meth_tot[i].first;
        N_hi = cumu_right_tot[i];
        k_hi = cumu_right_meth[i];
      }
      if (N_hi > 0 && N_low > 0) {
        double score = 0;
        const double p_hi = static_cast<double>(k_hi) / N_hi;
        const double p_low = static_cast<double>(k_low) / N_low;

        for (std::size_t j = 0; j < std::size(fg_alpha); ++j) {
          score += beta_max_likelihood(fg_alpha[j], fg_beta[j], bg_alpha[j],
                                       bg_beta[j], p_low, p_hi);
        }  // beta max likelihood using learned emissions
        score /= std::size(fg_alpha);
        if (p_hi > p_low && score > best_score) {
          best_idx = i;
          best_score = score;
        }
      }
    }
  return (best_score > -std::numeric_limits<double>::max())
           ? positions[best_idx]
           : std::numeric_limits<std::size_t>::max();
}

static void
get_boundary_positions(std::vector<Interval6> &bounds,
                       const std::vector<Interval6> &pmds,
                       const std::size_t &bin_size) {
  for (std::size_t i = 0; i < std::size(pmds); ++i) {
    bounds.push_back(pmds[i]);
    bounds.back().start =
      pmds[i].start > bin_size ? pmds[i].start - bin_size : 0;
    bounds.back().stop = pmds[i].start + bin_size;

    bounds.push_back(pmds[i]);
    bounds.back().start = pmds[i].stop > bin_size ? pmds[i].stop - bin_size : 0;
    bounds.back().stop = pmds[i].stop + bin_size;
  }
}

static void
get_optimized_boundary_likelihoods(
  const std::vector<std::string> &cpgs_file,
  const std::vector<Interval6> &bounds, const std::vector<bool> &array_status,
  const std::vector<double> &fg_alpha, const std::vector<double> &fg_beta,
  const std::vector<double> &bg_alpha, const std::vector<double> &bg_beta,
  std::vector<double> &boundary_scores,
  std::vector<std::size_t> &boundary_certainties) {
  // For weighting array contribution to boundary observations
  static constexpr double array_coverage_constant = 10;

  std::vector<bamxx::bgzf_file> in;
  for (std::size_t i = 0; i < std::size(cpgs_file); ++i) {
    in.emplace_back(cpgs_file[i], "r");
    if (get_has_counts_header(cpgs_file[i]))
      skip_counts_header(in.back());
  }

  std::map<std::size_t, std::pair<std::size_t, std::size_t>> pos_meth_tot;
  std::size_t n_meth{};
  std::size_t n_reads{};
  std::size_t bound_idx{};
  for (; bound_idx < std::size(bounds); ++bound_idx) {  // for each boundary
    for (std::size_t i = 0; i < std::size(in); ++i) {
      // get totals for all CpGs overlapping that boundary
      MSite site;
      while (read_site(in[i], site) &&
             !succeeds(site.chrom, site.pos, bounds[bound_idx])) {
        if (array_status[i])
          site.n_reads = array_coverage_constant;

        // check if CpG is inside boundary
        if (!precedes(site.chrom, site.pos, bounds[bound_idx])) {
          if (array_status[i]) {
            if (site.meth != -1) {
              n_meth = site.n_meth();
              n_reads = site.n_reads;
            }
            else {
              n_meth = 0;
              n_reads = 0;
            }
          }
          else {
            n_meth = site.n_meth();
            n_reads = site.n_reads;
          }
          if (pos_meth_tot.find(site.pos) == std::cend(pos_meth_tot))
            pos_meth_tot[site.pos] = std::make_pair(n_meth, n_reads);
          else {  // add this file's contribution to the site's methylation
            pos_meth_tot[site.pos].first += n_meth;
            pos_meth_tot[site.pos].second += n_reads;
          }
        }
      }
    }

    // Get the boundary position
    std::size_t boundary_position =
      (bounds[bound_idx].start + bounds[bound_idx].stop) / 2;

    std::size_t N_low = 0, k_low = 0, N_hi = 0, k_hi = 0;
    for (const auto &p : pos_meth_tot) {
      if (p.first < boundary_position) {
        N_low += p.second.second;
        k_low += p.second.first;
      }
      else {
        N_hi += p.second.second;
        k_hi += p.second.first;
      }
    }

    double score = 0;
    const double p_hi = static_cast<double>(k_hi) / N_hi;
    const double p_low = static_cast<double>(k_low) / N_low;

    if (bound_idx % 2) {  // its a right boundary, p_low should go with fg
      for (std::size_t j = 0; j < std::size(fg_alpha); ++j)
        score += beta_max_likelihood(fg_alpha[j], fg_beta[j], bg_alpha[j],
                                     bg_beta[j], p_low, p_hi);
    }
    else {  // its a left boundary, p_low should go with bg
      for (std::size_t j = 0; j < std::size(fg_alpha); ++j)
        score += beta_max_likelihood(bg_alpha[j], bg_beta[j], fg_alpha[j],
                                     fg_beta[j], p_low, p_hi);
    }
    boundary_certainties.push_back(std::min(N_low, N_hi));
    score /= std::size(fg_alpha);
    boundary_scores.push_back(exp(score));
    pos_meth_tot.clear();
  }
}

static void
find_exact_boundaries(const std::vector<std::string> &cpgs_file,
                      const std::vector<Interval6> &bounds,
                      const std::vector<bool> &array_status,
                      const std::vector<double> &fg_alpha,
                      const std::vector<double> &fg_beta,
                      const std::vector<double> &bg_alpha,
                      const std::vector<double> &bg_beta,
                      std::vector<std::size_t> &bound_site) {
  // For weighting array contribution to boundary observations
  static constexpr double array_coverage_constant = 10;

  std::vector<bamxx::bgzf_file> in;
  for (std::size_t i = 0; i < std::size(cpgs_file); ++i) {
    in.emplace_back(cpgs_file[i], "r");
    if (get_has_counts_header(cpgs_file[i]))
      skip_counts_header(in[i]);
  }

  std::map<std::size_t, std::pair<std::size_t, std::size_t>> pos_meth_tot;
  std::size_t n_meth{};
  std::size_t n_reads{};
  std::size_t bound_idx{};
  for (; bound_idx < std::size(bounds); ++bound_idx) {  // for each boundary
    for (std::size_t i = 0; i < std::size(in); ++i) {
      // get totals for all CpGs overlapping that boundary
      MSite site;
      while (read_site(in[i], site) &&
             !succeeds(site.chrom, site.pos, bounds[bound_idx])) {
        if (array_status[i])
          site.pos = array_coverage_constant;

        // check if CpG is inside boundary
        if (!precedes(site.chrom, site.pos, bounds[bound_idx])) {
          if (array_status[i]) {
            if (site.meth != -1) {
              n_meth = site.n_meth();
              n_reads = site.n_reads;
            }
            else {
              n_meth = 0;
              n_reads = 0;
            }
          }
          else {
            n_meth = site.n_meth();
            n_reads = site.n_reads;
          }
          auto it = pos_meth_tot.find(site.pos);
          if (it == end(pos_meth_tot)) {  // does not exist in map
            pos_meth_tot.emplace(site.pos, std::make_pair(n_meth, n_reads));
          }
          else {  // add this file's contribution to the CpG's methylation
            pos_meth_tot[site.pos].first += site.n_meth();
            pos_meth_tot[site.pos].second += site.n_reads;
          }
        }
      }
    }
    bound_site.push_back(find_best_bound(bound_idx % 2, pos_meth_tot, fg_alpha,
                                         fg_beta, bg_alpha, bg_beta));
    pos_meth_tot.clear();
  }
}

static void
optimize_boundaries(
  const std::size_t bin_size, const std::vector<std::string> &cpgs_file,
  std::vector<Interval6> &pmds, const std::vector<bool> &array_status,
  const std::vector<double> &fg_alpha, const std::vector<double> &fg_beta,
  const std::vector<double> &bg_alpha, const std::vector<double> &bg_beta) {
  std::vector<Interval6> bounds;
  get_boundary_positions(bounds, pmds, bin_size);
  std::vector<std::size_t> bound_site;
  find_exact_boundaries(cpgs_file, bounds, array_status, fg_alpha, fg_beta,
                        bg_alpha, bg_beta, bound_site);

  // Now reset the starts and ends of PMDs
  for (std::size_t i = 0; i < std::size(pmds); ++i) {
    const std::size_t start_site = bound_site[2 * i];
    if (start_site != std::numeric_limits<std::size_t>::max())
      pmds[i].start = start_site;
    const std::size_t end_site = bound_site[2 * i + 1];
    if (end_site != std::numeric_limits<std::size_t>::max())
      pmds[i].stop = end_site + 1;
  }

  // Now merge the pmds that are too close
  std::vector<std::size_t> dists;
  get_adjacent_distances(pmds, dists);
  std::sort(std::begin(dists), std::end(dists));

  // Need to use some randomization method here to figure out the
  // merging distance
  std::vector<std::pair<std::size_t, std::size_t>> dist_hist;
  std::size_t first = 0;
  for (std::size_t i = 1; i < std::size(dists); ++i)
    if (dists[i] != dists[i - 1]) {
      dist_hist.push_back(std::make_pair(dists[i - 1], i - first));
      first = i;
    }
  merge_nearby_pmd(2 * bin_size, pmds);

  // Last, get the cpg sites within 1 bin of each boundary and compute
  // the likelihood to get a "score" on the boundary
  bounds.clear();  // need updated boundaries after merging nearby PMDs
  get_boundary_positions(bounds, pmds, bin_size);

  std::vector<double> boundary_scores;
  std::vector<std::size_t> boundary_certainties;
  get_optimized_boundary_likelihoods(cpgs_file, bounds, array_status, fg_alpha,
                                     fg_beta, bg_alpha, bg_beta,
                                     boundary_scores, boundary_certainties);

  // Add the boundary scores to the PMD names
  // clang-format off
  for (std::size_t i = 0; i < std::size(pmds); ++i)
    pmds[i].name = pmds[i].name + ":" +
      std::to_string(boundary_scores[2 * i]) + ":" +
      std::to_string(boundary_certainties[2 * i]) + ":" +
      std::to_string(boundary_scores[2 * i + 1]) + ":" +
      std::to_string(boundary_certainties[2 * i + 1]);
  // clang-format on
}

double
get_score_cutoff_for_fdr(const std::vector<double> &scores, const double fdr) {
  if (fdr <= 0)
    return std::numeric_limits<double>::max();
  if (fdr > 1)
    return std::numeric_limits<double>::min();
  std::vector<double> local(scores);
  std::sort(std::begin(local), std::end(local));
  std::size_t i = 0;
  for (; i + 1 < std::size(local) &&
         local[i + 1] < fdr * static_cast<double>(i + 1) / std::size(local);
       ++i)
    ;
  return local[i] + 1.0 / std::size(scores);
}

static inline double
score_contribution(const std::pair<double, double> &m) {
  const double denom = m.first + m.second;
  return denom > 0 ? 1.0 - m.first / denom : 0.0;
}

static void
get_domain_scores(
  const std::vector<bool> &classes,
  const std::vector<std::vector<std::pair<double, double>>> &meth,
  const std::vector<std::size_t> &reset_points, std::vector<double> &scores) {
  const std::size_t n_replicates = std::size(meth);
  std::size_t reset_idx = 1;
  bool in_domain = false;
  double score = 0;
  for (std::size_t i = 0; i < std::size(classes); ++i) {
    if (reset_points[reset_idx] == i) {
      if (in_domain) {
        scores.push_back(score);
        score = 0;
        in_domain = false;
      }
      ++reset_idx;
    }
    if (classes[i]) {
      for (std::size_t r = 0; r < n_replicates; ++r)
        score += score_contribution(meth[r][i]);
      in_domain = true;
    }
    else if (in_domain) {
      scores.push_back(score);
      score = 0;
      in_domain = false;
    }
  }
  if (in_domain)
    scores.push_back(score);
}

static void
build_domains(const std::vector<Interval> &bins,
              const std::vector<std::size_t> &reset_points,
              const std::vector<bool> &classes,
              std::vector<Interval6> &domains) {
  std::size_t n_bins = 0, reset_idx = 1, prev_end = 0;
  bool in_domain = false;
  for (std::size_t i = 0; i < std::size(classes); ++i) {
    if (reset_points[reset_idx] == i) {
      if (in_domain) {
        domains.back().stop = prev_end;
        domains.back().score = n_bins;
        n_bins = 0;
        in_domain = false;
      }
      ++reset_idx;
    }
    if (classes[i]) {
      if (!in_domain) {
        domains.emplace_back(bins[i].chrom, bins[i].start, bins[i].stop,
                             std::string{}, 0, '+');
        in_domain = true;
      }
      ++n_bins;
    }
    else if (in_domain) {
      domains.back().stop = prev_end;
      domains.back().score = n_bins;
      n_bins = 0;
      in_domain = false;
    }
    prev_end = bins[i].stop;
  }
}

// Modified to take multiple replicates
template <class T, class U>
static void
separate_regions(const std::size_t desert_size,
                 std::vector<std::vector<Interval>> &bins,
                 std::vector<std::vector<T>> &meth,
                 std::vector<std::vector<U>> &reads,
                 std::vector<std::size_t> &reset_points,
                 std::vector<std::size_t> &dists_btwn_bins) {
  const std::size_t n_replicates = std::size(bins);

  // eliminate the zero-read cpg sites if no coverage in any rep
  std::size_t end_coord_of_prev = 0;
  std::size_t j = 0;
  for (std::size_t i = 0; i < std::size(bins[0]); ++i) {
    bool all_empty = true;
    std::size_t rep_idx = 0;
    while (all_empty && rep_idx < n_replicates) {
      if (reads[rep_idx][i] == 0)
        ++rep_idx;
      else
        all_empty = false;
    }
    if (!all_empty) {
      dists_btwn_bins.push_back(bins[0][i].start - end_coord_of_prev);
      end_coord_of_prev = bins[0][i].stop;
      for (std::size_t r = 0; r < n_replicates; ++r) {
        bins[r][j] = bins[r][i];
        meth[r][j] = meth[r][i];
        reads[r][j] = reads[r][i];
      }
      ++j;
    }
  }

  for (std::size_t r = 0; r < n_replicates; ++r) {
    bins[r].resize(j);
    meth[r].resize(j);
    reads[r].resize(j);
  }

  // segregate bins
  std::size_t prev_cpg = 0;
  for (std::size_t i = 0; i < std::size(bins[0]); ++i) {
    const std::size_t dist = (i > 0 && bins[0][i].chrom == bins[0][i - 1].chrom)
                               ? bins[0][i].start - prev_cpg
                               : std::numeric_limits<std::size_t>::max();
    if (dist > desert_size)
      reset_points.push_back(i);
    prev_cpg = bins[0][i].start;
  }
  assert(std::size(reset_points) > 0);
  reset_points.push_back(std::size(bins[0]));
}

static void
shuffle_bins(const std::size_t rng_seed, const TwoStateHMM &hmm,
             std::vector<std::vector<std::pair<double, double>>>
               meth,  // cppcheck-suppress[passedByValue]
             const std::vector<std::size_t> &reset_points,
             const std::vector<double> &start_trans,
             const std::vector<std::vector<double>> &trans,
             const std::vector<double> &end_trans,
             const std::vector<double> &fg_alpha,
             const std::vector<double> &fg_beta,
             const std::vector<double> &bg_alpha,
             const std::vector<double> &bg_beta,
             std::vector<double> &domain_scores,
             const std::vector<bool> &array_status) {
  const auto n_replicates = std::size(meth);
  auto eng = std::default_random_engine(rng_seed);
  for (std::size_t r = 0; r < n_replicates; ++r)
    std::shuffle(std::begin(meth[r]), std::end(meth[r]), eng);

  std::vector<bool> classes;
  std::vector<double> scores;
  hmm.PosteriorDecoding_rep(meth, reset_points, start_trans, trans, end_trans,
                            fg_alpha, fg_beta, bg_alpha, bg_beta, classes,
                            scores, array_status);
  get_domain_scores(classes, meth, reset_points, domain_scores);
  std::sort(std::begin(domain_scores), std::end(domain_scores));
}

static void
assign_p_values(const std::vector<double> &random_scores,
                const std::vector<double> &observed_scores,
                std::vector<double> &p_values) {
  const double n_randoms = std::max(std::size(random_scores), 1ul);
  for (auto scr : observed_scores) {
    const auto scr_itr =
      upper_bound(std::cbegin(random_scores), std::cend(random_scores), scr);
    p_values.push_back(std::distance(scr_itr, std::cend(random_scores)) /
                       n_randoms);
  }
}

static void
read_params_file(const bool verbose, const std::string &params_file,
                 double &fg_alpha, double &fg_beta, double &bg_alpha,
                 double &bg_beta, std::vector<double> &start_trans,
                 std::vector<std::vector<double>> &trans,
                 std::vector<double> &end_trans, double &fdr_cutoff) {
  std::string jnk;
  std::ifstream in(params_file.c_str());

  in >> jnk >> fg_alpha >> jnk >> fg_beta >> jnk >> bg_alpha >> jnk >>
    bg_beta >> jnk >> start_trans[0] >> jnk >> start_trans[1] >> jnk >>
    trans[0][0] >> jnk >> trans[0][1] >> jnk >> trans[1][0] >> jnk >>
    trans[1][1] >> jnk >> end_trans[0] >> jnk >> end_trans[1] >> jnk >>
    fdr_cutoff;

  if (verbose)
    std::cerr << "Read in params from " << params_file << '\n'
              << "FG_ALPHA\t" << fg_alpha << '\n'
              << "FG_BETA\t" << fg_beta << '\n'
              << "BG_ALPHA\t" << bg_alpha << '\n'
              << "BG_BETA\t" << bg_beta << '\n'
              << "S_F\t" << start_trans[0] << '\n'
              << "S_B\t" << start_trans[1] << '\n'
              << "F_F\t" << trans[0][0] << '\n'
              << "F_B\t" << trans[0][1] << '\n'
              << "B_F\t" << trans[1][0] << '\n'
              << "B_B\t" << trans[1][1] << '\n'
              << "F_E\t" << end_trans[0] << '\n'
              << "B_E\t" << end_trans[1] << '\n'
              << "FDR_CUTOFF\t" << fdr_cutoff << '\n';
}

static void
write_posteriors_file(const std::string &posteriors_file,
                      const std::vector<std::vector<Interval>> &bins,
                      const std::vector<double> &scores) {
  static constexpr auto decimal_precision = 10;
  std::ofstream out(posteriors_file);
  if (!out)
    throw std::runtime_error("failed to open: " + posteriors_file);
  out.precision(decimal_precision);
  for (std::size_t r = 0; r < std::size(scores); ++r)
    out << bins[0][r] << '\t' << scores[r] << '\n';
}

static void
write_params_file(const std::string &outfile,
                  const std::vector<double> &fg_alpha,
                  const std::vector<double> &fg_beta,
                  const std::vector<double> &bg_alpha,
                  const std::vector<double> &bg_beta,
                  const std::vector<double> &start_trans,
                  const std::vector<std::vector<double>> &trans,
                  const std::vector<double> &end_trans) {
  static constexpr auto decimal_precision = 30;
  std::ofstream out(outfile);
  if (!out)
    throw std::runtime_error("failed to open: " + outfile);
  out.precision(decimal_precision);
  // NOLINTBEGIN(*-avoid-magic-numbers)
  for (std::size_t r = 0; r < std::size(fg_alpha); ++r)
    out << "FG_ALPHA_" << r + 1 << "\t" << std::setw(14) << fg_alpha[r] << "\t"
        << "FG_BETA_" << r + 1 << "\t" << std::setw(14) << fg_beta[r] << "\t"
        << "BG_ALPHA_" << r + 1 << "\t" << std::setw(14) << bg_alpha[r] << "\t"
        << "BG_BETA_" << r + 1 << "\t" << std::setw(14) << bg_beta[r] << '\n';
  // NOLINTEND(*-avoid-magic-numbers)
  out << "S_F\t" << start_trans[0] << '\n'
      << "S_B\t" << start_trans[1] << '\n'
      << "F_F\t" << trans[0][0] << '\n'
      << "F_B\t" << trans[0][1] << '\n'
      << "B_F\t" << trans[1][0] << '\n'
      << "B_B\t" << trans[1][1] << '\n'
      << "F_E\t" << end_trans[0] << '\n'
      << "B_E\t" << end_trans[1] << '\n';
}

static bool
check_if_array_data(const std::string &infile) {
  bamxx::bgzf_file in(infile, "r");
  if (!in)
    throw std::runtime_error("bad file: " + infile);

  if (get_has_counts_header(infile))
    skip_counts_header(in);

  std::string line;
  getline(in, line);
  std::istringstream iss(line);
  std::string chrom, pos, strand, seq, meth, cov;
  iss >> chrom >> pos >> strand >> seq >> meth;
  return (!(iss >> cov));
}

static void
load_array_data(const std::size_t bin_size, const std::string &cpgs_file,
                std::vector<Interval> &bins,
                std::vector<std::pair<double, double>> &meth,
                std::vector<std::size_t> &reads) {
  // MAGIC. GS: minimum value for array?
  static constexpr double meth_min = 1.0e-2;

  bamxx::bgzf_file in(cpgs_file, "r");
  if (!in)
    throw std::runtime_error("bad sites file: " + cpgs_file);

  if (get_has_counts_header(cpgs_file))
    skip_counts_header(in);

  std::string curr_chrom;
  std::size_t prev_pos = 0ul, curr_pos = 0ul;
  double array_meth_bin = 0.0;
  double num_probes_in_bin = 0.0;

  MSite site;
  while (read_site(in, site)) {
    // TODO(MN): I think that the block below should be placed later in this
    // scope. At this location, the methylation level of the first site in a
    // new chrom is contributed to the last bin of the previous chrom.
    if (site.n_reads > 0) {  // its covered by a probe
      ++num_probes_in_bin;
      if (site.meth < meth_min)
        array_meth_bin += meth_min;
      else if (site.meth > 1.0 - meth_min)
        array_meth_bin += (1.0 - meth_min);
      else
        array_meth_bin += site.meth;
    }

    if (curr_chrom != site.chrom) {
      if (!curr_chrom.empty()) {
        if (site.chrom < curr_chrom)
          throw std::runtime_error("CpGs not sorted in file \"" + cpgs_file +
                                   "\"");
        bins.push_back(Interval(curr_chrom, curr_pos, prev_pos + 1));
        meth.push_back(std::make_pair(array_meth_bin, num_probes_in_bin));
        if (num_probes_in_bin > 0)
          reads.push_back(1);
        else
          reads.push_back(0);
      }
      curr_chrom = site.chrom;
      curr_pos = site.pos;
      array_meth_bin = 0.0;
      num_probes_in_bin = 0.0;
    }
    else if (curr_pos > site.pos) {
      throw std::runtime_error("CpGs not sorted in file \"" + cpgs_file + "\"");
    }
    else if (site.pos > curr_pos + bin_size) {
      bins.push_back(Interval(curr_chrom, curr_pos, curr_pos + bin_size));
      meth.push_back(std::make_pair(array_meth_bin, num_probes_in_bin));
      (num_probes_in_bin > 0) ? reads.push_back(1) : reads.push_back(0);

      array_meth_bin = 0.0;
      num_probes_in_bin = 0.0;
      curr_pos += bin_size;
      while (curr_pos + bin_size < site.pos) {
        bins.push_back(Interval(curr_chrom, curr_pos, curr_pos + bin_size));
        reads.push_back(0);
        meth.push_back(std::make_pair(0.0, 0.0));
        curr_pos += bin_size;
      }
    }
  }

  if (site.meth != -1) {  // its covered by a probe
    ++num_probes_in_bin;
    if (site.meth < meth_min)
      array_meth_bin += meth_min;
    else if (site.meth > 1.0 - meth_min)
      array_meth_bin += (1.0 - meth_min);
    else
      array_meth_bin += site.meth;
  }
  prev_pos = site.pos;
  if (!curr_chrom.empty()) {
    bins.push_back(Interval(curr_chrom, curr_pos, prev_pos + 1));
    meth.push_back(std::make_pair(array_meth_bin, num_probes_in_bin));
    if (num_probes_in_bin > 0)
      reads.push_back(1);
    else
      reads.push_back(0);
  }
}

static void
load_wgbs_data(const std::size_t bin_size, const std::string &cpgs_file,
               std::vector<Interval> &bins,
               std::vector<std::pair<double, double>> &meth,
               std::vector<std::size_t> &reads) {
  reads.clear();  // for safety
  meth.clear();
  bins.clear();

  // ADS: loading data each iteration should be put outside loop
  bamxx::bgzf_file in(cpgs_file, "r");
  if (!in)
    throw std::runtime_error("bad sites file: " + cpgs_file);

  if (get_has_counts_header(cpgs_file))
    skip_counts_header(in);

  // keep track of the chroms we've seen
  std::string curr_chrom;
  std::unordered_set<std::string> chroms_seen;

  MSite site;
  std::size_t prev_pos = 0ul;
  std::size_t sites_in_bin = 0ul;

  while (read_site(in, site)) {
    if (curr_chrom != site.chrom) {  // handle change of chrom
      if (sites_in_bin > 0)
        bins.back().stop = prev_pos;
      if (chroms_seen.find(site.chrom) != std::cend(chroms_seen))
        throw std::runtime_error("sites not sorted");
      chroms_seen.insert(site.chrom);
      curr_chrom = site.chrom;
      reads.push_back(0);
      meth.push_back(std::make_pair(0.0, 0.0));
      bins.push_back(Interval(site.chrom, 0, bin_size));
      sites_in_bin = 0;
    }
    prev_pos = site.pos;
    if (site.pos < bins.back().start)
      throw std::runtime_error("sites not sorted");
    while (bins.back().stop < site.pos) {
      sites_in_bin = 0;
      reads.push_back(0);
      meth.push_back(std::make_pair(0.0, 0.0));
      bins.push_back(
        Interval(site.chrom, bins.back().stop, bins.back().stop + bin_size));
    }
    reads.back() += site.n_reads;
    meth.back().first += site.n_meth();
    meth.back().second += site.n_unmeth();
    sites_in_bin++;
  }
  if (sites_in_bin > 0)
    bins.back().stop = prev_pos;
}

static void
remove_empty_bins_at_chrom_start(std::vector<Interval> &bins,
                                 std::vector<std::pair<double, double>> &meth,
                                 std::vector<std::size_t> &reads) {
  bool chrom_start = true;
  std::size_t j = 0;
  std::string prev_chrom = "";
  for (std::size_t i = 0; i < std::size(bins); ++i) {
    if (bins[i].chrom != prev_chrom) {
      chrom_start = true;
      prev_chrom = bins[i].chrom;
    }
    if (reads[i] > 0)
      chrom_start = false;
    if (!chrom_start) {
      reads[j] = reads[i];
      meth[j] = meth[i];
      bins[j] = bins[i];
      ++j;
    }
  }
  bins.erase(std::begin(bins) + j, std::end(bins));
  meth.erase(std::begin(meth) + j, std::end(meth));
  reads.erase(std::begin(reads) + j, std::end(reads));
}

static void
load_read_counts(const std::string &cpgs_file, const std::size_t bin_size,
                 std::vector<std::size_t> &reads) {
  reads.clear();  // for safety

  // ADS: loading data each iteration should be put outside loop
  bamxx::bgzf_file in(cpgs_file, "r");
  if (!in)
    throw std::runtime_error("bad methcounts file: " + cpgs_file);

  if (get_has_counts_header(cpgs_file))
    skip_counts_header(in);

  // keep track of where we are and what we've seen
  std::size_t bin_start = 0ul;
  std::string curr_chrom;
  std::unordered_set<std::string> chroms_seen;

  MSite site;
  while (read_site(in, site)) {
    if (curr_chrom != site.chrom) {  // handle change of chrom
      if (chroms_seen.find(site.chrom) != std::cend(chroms_seen))
        throw std::runtime_error("sites not sorted");
      chroms_seen.insert(site.chrom);
      bin_start = 0;
      curr_chrom = site.chrom;
      reads.push_back(0);
    }
    if (site.pos < bin_start)
      throw std::runtime_error("sites not sorted");
    for (; bin_start + bin_size < site.pos; bin_start += bin_size)
      reads.push_back(0);
    reads.back() += site.n_reads;
  }
}

static double
good_bins_frac(const std::vector<std::size_t> &cumulative,
               const std::size_t min_bin_size, const std::size_t bin_size,
               const std::size_t min_cov_to_pass) {
  // make sure the target bin size is a multiple of the minimum so we
  // have the resolution to construct the new bins
  assert(bin_size % min_bin_size == 0);

  // the step size corresponds to the number of minium sized bins that
  // would make up a new bin of the target size
  const std::size_t step_size = bin_size / min_bin_size;

  std::size_t passing_bins = 0, covered_bins = 0;

  std::size_t prev_total = 0;
  for (std::size_t i = 0; i + step_size < std::size(cumulative);
       i += step_size) {
    const std::size_t curr_cumulative = cumulative[i + step_size];
    const std::size_t bin_count = curr_cumulative - prev_total;
    covered_bins += (bin_count > 0);
    passing_bins += (bin_count >= min_cov_to_pass);
    prev_total = curr_cumulative;
  }

  // check if there is a leftover bin at the end
  if (std::size(cumulative) % step_size != 0) {
    const std::size_t bin_count = cumulative.back() - prev_total;
    covered_bins += (bin_count > 0);
    passing_bins += (bin_count >= min_cov_to_pass);
  }

  return static_cast<double>(passing_bins) / std::max(1ul, covered_bins);
}

static std::size_t
get_min_reads_for_confidence(const double conf_level) {
  // ADS: value of 0.5 below important; this is where the CI is widest
  static const double fixed_phat = 0.5;
  std::size_t n_reads = 0;
  double lower = 0.0, upper = 1.0;
  // ADS: should be doubling first, followed by bisection
  while (1.0 - conf_level < upper - lower) {
    ++n_reads;
    wilson_ci_for_binomial(1.0 - conf_level, n_reads, fixed_phat, lower, upper);
  }
  return n_reads;
}

// ADS: this function will return std::numeric_limits<std::size_t>::max() if the
// fraction of "good" bins is zero for all attempted bin sizes.
static std::size_t
binsize_selection(const std::size_t resolution, const std::size_t min_bin_sz,
                  const std::size_t max_bin_sz, const double conf_level,
                  const double min_frac_passed, const std::string &cpgs_file) {
  const std::size_t min_cov_to_pass = get_min_reads_for_confidence(conf_level);

  std::vector<std::size_t> reads;
  load_read_counts(cpgs_file, resolution, reads);

  std::partial_sum(std::cbegin(reads), std::cend(reads), std::begin(reads));

  double frac_passed = 0.0;
  std::size_t bin_size = min_bin_sz;

  while (bin_size < max_bin_sz && frac_passed < min_frac_passed) {
    frac_passed = good_bins_frac(reads, resolution, bin_size, min_cov_to_pass);
    if (frac_passed < min_frac_passed)
      bin_size += resolution;
  }
  return frac_passed < min_frac_passed ? std::numeric_limits<std::size_t>::max()
                                       : bin_size;
}

static void
load_bins(const std::size_t bin_size, const std::string &cpgs_file,
          std::vector<Interval> &bins,
          std::vector<std::pair<double, double>> &meth,
          std::vector<std::size_t> &reads, std::vector<bool> &array_status) {
  const bool is_array_data = check_if_array_data(cpgs_file);

  array_status.push_back(is_array_data);

  if (is_array_data)
    load_array_data(bin_size, cpgs_file, bins, meth, reads);
  else {
    load_wgbs_data(bin_size, cpgs_file, bins, meth, reads);
    remove_empty_bins_at_chrom_start(bins, meth, reads);
  }
}

static void
get_union_of_bins(const std::vector<std::vector<Interval>> &orig,
                  std::vector<Interval> &bins) {
  const auto overlaps = [](const auto &x, const auto &y) {
    const auto cmp = x.chrom.compare(y.chrom);
    return cmp == 0 && std::max(x.start, x.start) < std::min(x.stop, y.stop);
  };

  // flatten the set of sorted bins
  bins.clear();
  for (auto &&i : orig)
    bins.insert(std::end(bins), std::cbegin(i), std::cend(i));

  // merge each sorted interval of bins
  const auto first = std::begin(bins);
  auto middle = std::begin(bins);
  for (std::size_t i = 1; i < std::size(orig); ++i) {
    middle += std::size(orig[i - 1]);
    std::inplace_merge(first, middle, middle + std::size(orig[i]));
  }
  // ensure unique bins
  bins.erase(std::unique(std::begin(bins), std::end(bins)), std::end(bins));
  bins.shrink_to_fit();

  // make sure all bins are aligned at same boundaries
  for (std::size_t i = 1; i < std::size(bins); ++i)
    if (overlaps(bins[i - 1], bins[i]))
      throw std::runtime_error("bins from reps not aligned");
}

static void
add_missing_bins(const std::vector<Interval> &all_bins,
                 std::vector<Interval> &bins,
                 std::vector<std::pair<double, double>> &meth) {
  const std::size_t n_bins = std::size(all_bins);
  std::vector<std::pair<double, double>> tmp_meth(n_bins);

  std::size_t j = 0;  // assume j range no larger than i range
  for (std::size_t i = 0; i < n_bins; ++i) {
    if (all_bins[i] == bins[j])
      tmp_meth[i] = meth[j++];
    else
      tmp_meth[i] = {0.0, 0.0};
  }
  std::swap(meth, tmp_meth);
  bins = all_bins;
}

static void
write_empty_summary(const std::string &summary_file) {
  if (summary_file.empty())
    return;
  std::ofstream summary_out(summary_file);
  if (!summary_out)
    throw std::runtime_error("failed to open: " + summary_file);
  summary_out << to_string(pmd_summary({})) << '\n';
}

int
main_pmd(int argc, char *argv[]) {  // NOLINT(*-avoid-c-arrays)
  try {
    static constexpr auto min_observations_for_inference{100ul};
    static constexpr auto default_self_transition{0.99};
    static constexpr auto default_transition{0.01};
    static constexpr auto default_start_transition{0.5};
    static constexpr auto default_end_transition{1e-10};

    // magic numbers from paper: highest jaccard index to wgbs
    static constexpr auto bin_size_for_array{1000ul};
    static constexpr auto desert_size_for_array{200000ul};

    // MAGIC: corrections for small values (not parameters):
    static constexpr auto tolerance{1e-5};
    static constexpr auto min_prob{1e-10};

    static constexpr auto default_fg_alpha{0.05};
    static constexpr auto default_fg_beta{0.95};
    static constexpr auto default_bg_alpha{0.95};
    static constexpr auto default_bg_beta{0.05};
    const auto init_trans = [&] {
      return std::vector<std::vector<double>>{
        {default_self_transition, default_transition},
        {default_transition, default_self_transition},
      };
    };

    // NOLINTBEGIN(*-avoid-magic-numbers)
    std::size_t resolution = 500;
    std::size_t rng_seed = 408;
    std::size_t desert_size = 5000;
    std::size_t bin_size = 1000;
    std::size_t max_iterations = 10;
    // NOLINTEND(*-avoid-magic-numbers)

    // run mode flags
    bool DEBUG = false;
    bool verbose = false;
    bool ARRAY_MODE = false;
    bool fixed_bin_size = false;

    std::string summary_file;
    std::string params_in_files;
    std::string params_out_file;
    std::string posteriors_out_prefix;
    std::string outfile;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(argv[0],  // NOLINT(*-pointer-arithmetic)
                           description, "<methcount-files>");
    opt_parse.add_opt("out", 'o', "output file", true, outfile);
    opt_parse.add_opt("desert", 'd', "max dist between bins with data in PMD",
                      false, desert_size);
    opt_parse.add_opt("fixedbin", 'f', "Fixed bin size", false, fixed_bin_size);
    opt_parse.add_opt("bin", 'b', "Starting bin size", false, bin_size);
    opt_parse.add_opt("arraymode", 'a', "All samples are array", false,
                      ARRAY_MODE);
    opt_parse.add_opt("itr", 'i', "max iterations", false, max_iterations);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, verbose);
    opt_parse.add_opt("debug", 'D', "print more run info", false, DEBUG);
    opt_parse.add_opt("params-in", 'P',
                      "HMM parameter files for "
                      "individual methylomes (separated with comma)",
                      false, params_in_files);
    opt_parse.add_opt("posteriors-out", 'r',
                      "write out posterior probabilities in methcounts format",
                      false, posteriors_out_prefix);
    opt_parse.add_opt("summary", 'S', "write summary output here", false,
                      summary_file);
    opt_parse.add_opt("params-out", 'p', "write HMM parameters to this file",
                      false, params_out_file);
    opt_parse.add_opt("seed", 's', "specify random seed", false, rng_seed);

    std::vector<std::string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      std::cerr << opt_parse.help_message() << '\n'
                << opt_parse.about_message() << '\n';
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      std::cerr << opt_parse.about_message() << '\n';
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      std::cerr << opt_parse.option_missing_message() << '\n';
      return EXIT_SUCCESS;
    }
    if (leftover_args.empty()) {
      std::cerr << opt_parse.help_message() << '\n';
      return EXIT_SUCCESS;
    }
    /****************** END COMMAND LINE OPTIONS *****************/

    resolution = std::min(bin_size, resolution);

    std::vector<std::string> cpgs_file = leftover_args;
    std::vector<std::string> params_in_file;
    if (!params_in_files.empty())
      params_in_file = smithlab::split(params_in_files, ",", false);

    const std::size_t n_replicates = std::size(cpgs_file);
    std::for_each(std::cbegin(cpgs_file), std::cend(cpgs_file),
                  [](const auto &x) {
                    if (!is_msite_file(x))
                      throw std::runtime_error("malformed counts file: " + x);
                  });

    bool insufficient_data = false;  // ADS: this is used now to detect
                                     // when the counts files have
                                     // lines for CpG sites, but no
                                     // counts.

    // Sanity checks input file format and dynamically selects bin
    // size from WGBS samples.
    if (!fixed_bin_size && !ARRAY_MODE) {
      static constexpr auto max_bin_size{500000ul};
      static constexpr auto min_bin_size{1000ul};
      static constexpr auto desert_size_multiplier{5};  // ADS: what is this?
      static constexpr auto default_confidence_interval{0.80};
      static constexpr auto default_prop_accept{0.80};
      if (verbose)
        std::cerr << "[DYNAMICALLY SELECTING BIN SIZE]" << '\n';
      double confidence_interval = default_confidence_interval;
      double prop_accept = default_prop_accept;
      for (std::size_t i = 0; i < n_replicates && !insufficient_data; ++i) {
        const bool arrayData = check_if_array_data(cpgs_file[i]);
        if (!arrayData) {
          bin_size =
            binsize_selection(resolution, min_bin_size, max_bin_size,
                              confidence_interval, prop_accept, cpgs_file[i]);
          if (bin_size == std::numeric_limits<std::size_t>::max())
            insufficient_data = true;
          desert_size = desert_size_multiplier * bin_size;
        }
        else {
          bin_size = bin_size_for_array;
          desert_size = desert_size_for_array;
        }
      }
    }
    else if (ARRAY_MODE) {
      bin_size = bin_size_for_array;
      desert_size = desert_size_for_array;
    }
    else {
      desert_size = std::max(desert_size, bin_size);
    }

    if (insufficient_data) {
      // ADS: first check for insufficient data; another is needed if
      // fixed bin size is used
      if (verbose)
        std::cerr << "EXITING: INSUFFICIENT DATA" << '\n';
      if (!summary_file.empty())
        write_empty_summary(summary_file);
      return EXIT_SUCCESS;
    }

    if (verbose)
      std::cerr << "[READING IN AT BIN SIZE " << bin_size << "]\n";

    // separate the regions by chrom and by desert
    std::vector<std::vector<Interval>> bins(n_replicates);
    std::vector<std::vector<std::pair<double, double>>> meth(n_replicates);
    std::vector<std::vector<std::size_t>> reads(n_replicates);
    std::vector<bool> array_status;

    for (std::size_t i = 0; i < n_replicates && !insufficient_data; ++i) {
      if (verbose)
        std::cerr << "[READING CPGS AND METH PROPS] from " << cpgs_file[i]
                  << '\n';

      load_bins(bin_size, cpgs_file[i], bins[i], meth[i], reads[i],
                array_status);
      const double total_observations =
        std::accumulate(std::cbegin(reads[i]), std::cend(reads[i]), 0.0);
      if (total_observations <= std::numeric_limits<double>::min())
        insufficient_data = true;
      if (verbose)
        std::cerr << "TOTAL BINS: " << std::size(bins[i]) << '\n'
                  << "MEAN COVERAGE: "
                  << total_observations / std::max(1ul, std::size(reads[i]))
                  << '\n';
    }

    if (insufficient_data) {
      // ADS: second check for insufficient data; another is needed if
      // filtered number of bins is too few
      if (verbose)
        std::cerr << "EXITING: INSUFFICIENT DATA" << '\n';
      if (!summary_file.empty())
        write_empty_summary(summary_file);
      return EXIT_SUCCESS;
    }

    if (n_replicates > 1) {
      bool need_to_adjust_bins = false;
      const std::size_t n_bins_0 = std::size(bins[0]);
      for (std::size_t i = 1; i < n_replicates && !need_to_adjust_bins; ++i)
        if (std::size(bins[i]) != n_bins_0)
          need_to_adjust_bins = true;

      if (need_to_adjust_bins) {
        std::vector<Interval> all_bins;
        get_union_of_bins(bins, all_bins);
        for (std::size_t i = 0; i < std::size(bins); ++i)
          add_missing_bins(all_bins, bins[i], meth[i]);
      }
    }

    // separate regions by chrom and desert; eliminate isolated Bins
    std::vector<std::size_t> reset_points;
    std::vector<std::size_t> dists_btwn_bins;
    if (verbose)
      std::cerr << "[separating by CpG desert]" << '\n';
    separate_regions(desert_size, bins, meth, reads, reset_points,
                     dists_btwn_bins);
    if (size(bins[0]) < min_observations_for_inference)
      insufficient_data = true;

    if (insufficient_data) {
      // ADS: final check for sufficient data failed; too few bins
      // after filtering
      if (verbose)
        std::cerr << "EXITING: INSUFFICIENT DATA" << '\n';
      if (!summary_file.empty())
        write_empty_summary(summary_file);
      return EXIT_SUCCESS;
    }

    if (verbose)
      std::cerr << "bins retained: " << std::size(bins[0]) << '\n'
                << "number of distances between: " << std::size(dists_btwn_bins)
                << '\n'
                << "deserts removed: " << size(reset_points) - 2 << '\n';

    /****************** Read in params *****************/
    auto start_trans = std::vector<double>(2, default_start_transition);
    auto end_trans = std::vector<double>(2, default_end_transition);
    auto trans = init_trans();
    const TwoStateHMM hmm(min_prob, tolerance, max_iterations, verbose, DEBUG);
    std::vector<double> reps_fg_alpha(n_replicates, default_fg_alpha);
    std::vector<double> reps_fg_beta(n_replicates, default_fg_beta);
    std::vector<double> reps_bg_alpha(n_replicates, default_bg_alpha);
    std::vector<double> reps_bg_beta(n_replicates, default_bg_beta);
    double score_cutoff_for_fdr = std::numeric_limits<double>::max();
    if (!params_in_file.empty()) {
      // read parameters files
      for (auto i = 0u; i < n_replicates; ++i)
        read_params_file(verbose, params_in_file[i], reps_fg_alpha[i],
                         reps_fg_beta[i], reps_bg_alpha[i], reps_bg_beta[i],
                         start_trans, trans, end_trans, score_cutoff_for_fdr);
    }

    // train model (default behavior; not done when params supplied)
    if (max_iterations > 0)
      hmm.BaumWelchTraining_rep(meth, reset_points, start_trans, trans,
                                end_trans, reps_fg_alpha, reps_fg_beta,
                                reps_bg_alpha, reps_bg_beta, array_status);

    if (!params_out_file.empty()) {
      // write all the HMM parameters
      write_params_file(params_out_file, reps_fg_alpha, reps_fg_beta,
                        reps_bg_alpha, reps_bg_beta, start_trans, trans,
                        end_trans);
    }

    /***********************************/

    if (!posteriors_out_prefix.empty()) {
      std::vector<double> into_scores;
      hmm.TransitionPosteriors_rep(meth, reset_points, start_trans, trans,
                                   end_trans, reps_fg_alpha, reps_fg_beta,
                                   reps_bg_alpha, reps_bg_beta, array_status, 2,
                                   into_scores);
      write_posteriors_file(posteriors_out_prefix + ".intoTrans", bins,
                            into_scores);
      std::vector<double> outof_scores;
      hmm.TransitionPosteriors_rep(meth, reset_points, start_trans, trans,
                                   end_trans, reps_fg_alpha, reps_fg_beta,
                                   reps_bg_alpha, reps_bg_beta, array_status, 1,
                                   outof_scores);
      write_posteriors_file(posteriors_out_prefix + ".outofTrans", bins,
                            outof_scores);
    }

    std::vector<bool> classes;
    std::vector<double> scores;
    hmm.PosteriorDecoding_rep(meth, reset_points, start_trans, trans, end_trans,
                              reps_fg_alpha, reps_fg_beta, reps_bg_alpha,
                              reps_bg_beta, classes, scores, array_status);

    if (!posteriors_out_prefix.empty())
      write_posteriors_file(posteriors_out_prefix + ".posteriors", bins,
                            scores);

    std::vector<double> domain_scores;
    get_domain_scores(classes, meth, reset_points, domain_scores);

    if (verbose)
      std::cerr << "[RANDOMIZING SCORES FOR FDR]" << '\n';

    std::vector<double> random_scores;
    shuffle_bins(rng_seed, hmm, meth, reset_points, start_trans, trans,
                 end_trans, reps_fg_alpha, reps_fg_beta, reps_bg_alpha,
                 reps_bg_beta, random_scores, array_status);

    std::vector<double> p_values;
    assign_p_values(random_scores, domain_scores, p_values);

    if (score_cutoff_for_fdr == std::numeric_limits<double>::max() &&
        !p_values.empty()) {
      static constexpr auto default_p_value_cutoff{0.01};
      score_cutoff_for_fdr =
        get_score_cutoff_for_fdr(p_values, default_p_value_cutoff);
    }
    if (!params_out_file.empty()) {
      std::ofstream out(params_out_file, std::ios::app);
      if (!out)
        throw std::runtime_error("failed to open: " + params_out_file);
      out << "FDR_CUTOFF\t" << score_cutoff_for_fdr << '\n';
    }

    std::vector<Interval6> domains;
    build_domains(bins[0], reset_points, classes, domains);

    std::size_t good_pmd_count = 0;
    std::vector<Interval6> good_domains;
    for (std::size_t i = 0; i < std::size(domains); ++i)
      if (p_values[i] < score_cutoff_for_fdr) {
        good_domains.push_back(domains[i]);
        good_domains.back().name = "PMD" + std::to_string(good_pmd_count++);
      }

    optimize_boundaries(bin_size, cpgs_file, good_domains, array_status,
                        reps_fg_alpha, reps_fg_beta, reps_bg_alpha,
                        reps_bg_beta);

    if (!summary_file.empty()) {
      std::ofstream summary_out(summary_file);
      if (!summary_out)
        throw std::runtime_error("failed to open: " + summary_file);
      summary_out << to_string(pmd_summary(good_domains)) << '\n';
    }

    std::ofstream out(outfile);
    if (!out)
      throw std::runtime_error("failed to open: " + outfile);
    std::copy(std::cbegin(good_domains), std::cend(good_domains),
              std::ostream_iterator<Interval6>(out, "\n"));
  }
  catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

// NOLINTEND(*-narrowing-conversions)
