/* Copyright (C) 2018 University of Southern California
 *                         Andrew D Smith
 * Authors: Andrew D. Smith, Song Qiang, Jenny Qu,
 *          Benjamin Decato, Guilherme Sena
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

#include <gsl/gsl_sf.h>

#include <bamxx.hpp>

#include <numeric>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>
#include <stdexcept>
#include <unordered_set>
#include <random>
#include <sstream>

#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"
#include "OptionParser.hpp"
#include "bsutils.hpp"
#include "counts_header.hpp"

#include "TwoStateHMM_PMD.hpp"
#include "MSite.hpp"

using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::cerr;
using std::numeric_limits;
using std::max;
using std::min;
using std::make_pair;
using std::pair;
using std::runtime_error;

using std::ostream;
using std::ofstream;
using std::to_string;

using bamxx::bgzf_file;

template<typename T> using num_lim = std::numeric_limits<T>;

struct pmd_summary {
  pmd_summary(const vector<GenomicRegion> &pmds) {
    pmd_count = pmds.size();
    pmd_total_size = accumulate(cbegin(pmds), cend(pmds), 0ul,
                                [](const uint64_t t, const GenomicRegion &p) {
                                  return t + p.get_width(); });
    pmd_mean_size =
      static_cast<double>(pmd_total_size)/std::max(pmd_count, static_cast<uint64_t>(1));
  }
  // pmd_count is the number of identified PMDs.
  uint64_t pmd_count{};
  // total_pmd_size is the sum of the sizes of the identified PMDs
  uint64_t pmd_total_size{};
  // mean_pmd_size is the mean size of the identified PMDs
  double pmd_mean_size{};

  string tostring() {
    std::ostringstream oss;
    oss << "pmd_count: " << pmd_count << endl
        << "pmd_total_size: " << pmd_total_size << endl
        << "pmd_mean_size: "
        << std::fixed << std::setprecision(2) << pmd_mean_size;

    return oss.str();
  }
};

static void
get_adjacent_distances(const vector<GenomicRegion> &pmds,
                       vector<size_t> &dists) {
  for (size_t i = 1; i < pmds.size(); ++i)
    if (pmds[i].same_chrom(pmds[i - 1]))
      dists.push_back(pmds[i].get_start() - pmds[i - 1].get_end());
}


static bool
precedes(const string &chrom, const size_t position,
         const GenomicRegion &r) {
  return (chrom < r.get_chrom() ||
          (chrom == r.get_chrom() &&
           position < r.get_start()));
}


static bool
succeeds(const string &chrom, const size_t position,
         const GenomicRegion &r) {
  return (r.get_chrom() < chrom ||
          (chrom == r.get_chrom() &&
           r.get_end() <= position));
}


static void
merge_nearby_pmd(const size_t max_merge_dist,
                 vector<GenomicRegion> &pmds) {
  size_t j = 0;
  for (size_t i = 1; i < pmds.size(); ++i) {
    if (pmds[j].same_chrom(pmds[i]) &&
        pmds[i].get_start() - pmds[j].get_end() <= max_merge_dist) {
      pmds[j].set_end(pmds[i].get_end());
      const string combined_name(pmds[j].get_name() + pmds[i].get_name());
      pmds[j].set_name(combined_name);
    }
    else pmds[++j] = pmds[i];
  }
  pmds.resize(j);
}

inline double
beta_max_likelihood(const double fg_alpha, const double fg_beta,
                    const double bg_alpha, const double bg_beta,
                    const double p_low, const double p_hi) {
  return (fg_alpha - 1.0)*log(p_low) +
    (fg_beta - 1.0)*log(1.0 - p_low) -
    gsl_sf_lnbeta(fg_alpha, fg_beta) +
    (bg_alpha - 1.0)*log(p_hi) +
    (bg_beta - 1.0)*log(1.0 - p_hi) -
    gsl_sf_lnbeta(bg_alpha, bg_beta);
}

static size_t
find_best_bound(const bool IS_RIGHT_BOUNDARY,
                std::map<size_t, pair<size_t, size_t> > &pos_meth_tot,
                const vector<double> &fg_alpha,
                const vector<double> &fg_beta,
                const vector<double> &bg_alpha,
                const vector<double> &bg_beta) {

  vector<pair<size_t,size_t> > meth_tot;
  vector<size_t> positions;
  for (auto &&it : pos_meth_tot) {
    positions.push_back(it.first);
    meth_tot.push_back(it.second);
  }

  vector<size_t> cumu_left_meth(meth_tot.size(), 0);
  vector<size_t> cumu_left_tot(meth_tot.size(), 0);
  vector<size_t> cumu_right_meth(meth_tot.size(), 0);
  vector<size_t> cumu_right_tot(meth_tot.size(), 0);
  if (meth_tot.size() > 0)
    for (size_t i = 1; i < meth_tot.size()-1; ++i) {
      const size_t j = meth_tot.size() - 1 - i;
      cumu_left_meth[i] = cumu_left_meth[i-1] + meth_tot[i-1].first;
      cumu_left_tot[i] = cumu_left_tot[i-1] + meth_tot[i-1].second;
      cumu_right_meth[j] = cumu_right_meth[j+1] + meth_tot[j+1].first;
      cumu_right_tot[j] = cumu_right_tot[j+1] + meth_tot[j+1].second;
    }

  size_t best_idx = 0;
  double best_score = -num_lim<double>::max();
  if (meth_tot.size() > 0)
    for (size_t i = 1; i < meth_tot.size()-1; ++i) {
      size_t N_low, k_low, N_hi, k_hi;
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
        const double p_hi = static_cast<double>(k_hi)/N_hi;
        const double p_low = static_cast<double>(k_low)/N_low;

        for (size_t j = 0; j < fg_alpha.size(); ++j) {
          score += beta_max_likelihood(fg_alpha[j], fg_beta[j],
                                       bg_alpha[j], bg_beta[j],
                                       p_low, p_hi);
        } // beta max likelihood using learned emissions
        score /= fg_alpha.size();
        if (p_hi > p_low && score > best_score) {
          best_idx = i;
          best_score = score;
        }
      }
    }
  return (best_score > -num_lim<double>::max()) ?
    positions[best_idx] : num_lim<size_t>::max();
}


static void
get_boundary_positions(vector<GenomicRegion> &bounds,
                       const vector<GenomicRegion> &pmds,
                       const size_t &bin_size) {
  for (size_t i = 0; i < pmds.size(); ++i) {
    bounds.push_back(pmds[i]);
    pmds[i].get_start() > bin_size ?
      bounds.back().set_start(pmds[i].get_start() - bin_size) :
      bounds.back().set_start(0);
    bounds.back().set_end(pmds[i].get_start() + bin_size);

    bounds.push_back(pmds[i]);
    pmds[i].get_end() > bin_size ?
      bounds.back().set_start(pmds[i].get_end() - bin_size) :
      bounds.back().set_start(0);
    bounds.back().set_end(pmds[i].get_end() + bin_size);
  }
}

static void
get_optimized_boundary_likelihoods(const vector<string> &cpgs_file,
                                   vector<GenomicRegion> &bounds,
                                   const vector<bool> &array_status,
                                   const vector<double> &fg_alpha,
                                   const vector<double> &fg_beta,
                                   const vector<double> &bg_alpha,
                                   const vector<double> &bg_beta,
                                   vector<double> &boundary_scores,
                                   vector<size_t> &boundary_certainties) {
  // MAGIC NUMBER FOR WEIGHTING ARRAY
  // CONTRIBUTION TO BOUNDARY OBSERVATIONS
  static const double array_coverage_constant = 10;

  vector<bgzf_file*> in(cpgs_file.size());
  for (size_t i = 0; i < cpgs_file.size(); ++i) {
    in[i] = new bgzf_file(cpgs_file[i], "r");
    if (get_has_counts_header(cpgs_file[i]))
      skip_counts_header(*in[i]);
  }

  std::map<size_t, pair<size_t, size_t> > pos_meth_tot;
  size_t n_meth = 0ul, n_reads = 0ul;
  size_t bound_idx = 0;
  for (; bound_idx < bounds.size(); ++bound_idx) { // for each boundary
    for (size_t i = 0; i < in.size(); ++i) {
      // get totals for all CpGs overlapping that boundary

      MSite site;
      while (read_site(*in[i], site) &&
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
          auto it(pos_meth_tot.find(site.pos));
          if (it == end(pos_meth_tot))
            pos_meth_tot[site.pos] = make_pair(n_meth, n_reads);

          else { // add this file's contribution to the site's methylation
            pos_meth_tot[site.pos].first += n_meth;
            pos_meth_tot[site.pos].second += n_reads;
          }
        }
      }
    }

    // Get the boundary position
    size_t boundary_position =
      (bounds[bound_idx].get_start() + bounds[bound_idx].get_end())/2;

    size_t N_low = 0, k_low = 0, N_hi = 0, k_hi = 0;
    for (auto it = begin(pos_meth_tot); it != end(pos_meth_tot); ++it) {
      if (it->first < boundary_position) {
        N_low += it->second.second;
        k_low += it->second.first;
      }
      else{
        N_hi += it->second.second;
        k_hi += it->second.first;
      }
    }

    double score = 0;
    const double p_hi = static_cast<double>(k_hi)/N_hi;
    const double p_low = static_cast<double>(k_low)/N_low;

    if (bound_idx % 2) { // its a right boundary, p_low should go with fg
      for (size_t j = 0; j < fg_alpha.size(); ++j)
        score += beta_max_likelihood(fg_alpha[j], fg_beta[j],
                                     bg_alpha[j], bg_beta[j],
                                     p_low, p_hi);
    }
    else { // its a left boundary, p_low should go with bg
      for (size_t j = 0; j < fg_alpha.size(); ++j)
        score += beta_max_likelihood(bg_alpha[j], bg_beta[j],
                                     fg_alpha[j], fg_beta[j],
                                     p_low, p_hi);


    }
    boundary_certainties.push_back(std::min(N_low,N_hi));
    score /= fg_alpha.size();
    boundary_scores.push_back(exp(score));
    pos_meth_tot.clear();
  }

  for (auto &&fp : in)
    delete fp;
}


static void
find_exact_boundaries(const vector<string> &cpgs_file,
                      vector<GenomicRegion> &bounds,
                      const vector<bool> &array_status,
                      const vector<double> &fg_alpha,
                      const vector<double> &fg_beta,
                      const vector<double> &bg_alpha,
                      const vector<double> &bg_beta,
                      vector<size_t> &bound_site) {
  // MAGIC NUMBER FOR WEIGHTING ARRAY
  // CONTRIBUTION TO BOUNDARY OBSERVATIONS
  static const double array_coverage_constant = 10;

  vector<bgzf_file*> in(cpgs_file.size());
  for (size_t i = 0; i < cpgs_file.size(); ++i) {
    in[i] = new bgzf_file(cpgs_file[i], "r");
    if (get_has_counts_header(cpgs_file[i]))
      skip_counts_header(*in[i]);
  }

  std::map<size_t, pair<size_t, size_t> > pos_meth_tot;
  size_t n_meth = 0ul, n_reads = 0ul;
  size_t bound_idx = 0;
  for (; bound_idx < bounds.size(); ++bound_idx) { // for each boundary
    for (size_t i = 0; i < in.size(); ++i) {
      // get totals for all CpGs overlapping that boundary

      MSite site;
      while (read_site(*in[i], site) &&
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
          if (it == end(pos_meth_tot)) {// does not exist in map
            pos_meth_tot.emplace(site.pos, make_pair(n_meth, n_reads));
          }
          else { // add this file's contribution to the CpG's methylation
            pos_meth_tot[site.pos].first += site.n_meth();
            pos_meth_tot[site.pos].second += site.n_reads;
          }
        }
      }
    }
    bound_site.push_back(find_best_bound(bound_idx % 2, pos_meth_tot,
                                         fg_alpha, fg_beta,
                                         bg_alpha, bg_beta));
    pos_meth_tot.clear();
  }
  for (size_t i = 0; i < in.size(); ++i)
    delete in[i];
}


static void
optimize_boundaries(const size_t bin_size,
                    const vector<string> &cpgs_file,
                    vector<GenomicRegion> &pmds,
                    const vector<bool> &array_status,
                    const vector<double> &fg_alpha,
                    const vector<double> &fg_beta,
                    const vector<double> &bg_alpha,
                    const vector<double> &bg_beta) {

  vector<GenomicRegion> bounds;
  get_boundary_positions(bounds, pmds, bin_size);
  vector<size_t> bound_site;
  find_exact_boundaries(cpgs_file, bounds, array_status, fg_alpha,
                        fg_beta, bg_alpha, bg_beta,
                        bound_site);

  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  ///// NOW RESET THE STARTS AND ENDS OF PMDs
  /////
  for (size_t i = 0; i < pmds.size(); ++i) {
    const size_t start_site = bound_site[2*i];
    if (start_site != num_lim<size_t>::max())
      pmds[i].set_start(start_site);
    const size_t end_site = bound_site[2*i + 1];
    if (end_site != num_lim<size_t>::max())
      pmds[i].set_end(end_site + 1);
  }

  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  ///// NOW MERGE THE PMDS THAT ARE TOO CLOSE
  /////
  vector<size_t> dists;
  get_adjacent_distances(pmds, dists);
  sort(begin(dists), end(dists));

  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  // NEED TO USE SOME RANDOMIZATION METHOD HERE TO FIGURE OUT THE
  // MERGING DISTANCE
  vector<pair<size_t, size_t> > dist_hist;
  size_t first = 0;
  for (size_t i = 1; i < dists.size(); ++i)
    if (dists[i] != dists[i - 1]) {
      dist_hist.push_back(make_pair(dists[i - 1], i - first));
      first = i;
    }

  merge_nearby_pmd(2*bin_size, pmds);

  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  // LAST, GET THE CPG SITES WITHIN 1 BIN OF EACH BOUNDARY AND COMPUTE
  // THE LIKELIHOOD TO GET A "SCORE" ON THE BOUNDARY
  ///////
  bounds.clear(); // need updated boundaries after merging nearby PMDs
  get_boundary_positions(bounds, pmds, bin_size);

  vector<double> boundary_scores;
  vector<size_t> boundary_certainties;
  get_optimized_boundary_likelihoods(cpgs_file, bounds, array_status,
                                     fg_alpha, fg_beta, bg_alpha,
                                     bg_beta, boundary_scores,
                                     boundary_certainties);

  // Add the boundary scores to the PMD names
  for (size_t i = 0; i < pmds.size(); ++i)
    pmds[i].set_name(pmds[i].get_name()
                     + ":" + to_string(boundary_scores[2*i])
                     + ":" + to_string(boundary_certainties[2*i])
                     + ":" + to_string(boundary_scores[2*i+1])
                     + ":" + to_string(boundary_certainties[2*i+1]));
}


double
get_score_cutoff_for_fdr(const vector<double> &scores, const double fdr) {
  if (fdr <= 0)
    return num_lim<double>::max();
  else if (fdr > 1)
    return num_lim<double>::min();
  vector<double> local(scores);
  std::sort(begin(local), end(local));
  size_t i = 0;
  for (; i < local.size() - 1 &&
         local[i+1] < fdr*static_cast<double>(i+1)/local.size(); ++i);
  return local[i] + 1.0/scores.size();
}


static inline double
score_contribution(const pair<double, double> &m) {
  const double denom = m.first + m.second;
  return (denom > 0) ? 1.0 - m.first/denom : 0.0;
}


static void
get_domain_scores(const vector<bool> &classes,
                  const vector<vector<pair<double, double> > > &meth,
                  const vector<size_t> &reset_points,
                  vector<double> &scores) {

  const size_t n_replicates = meth.size();
  size_t reset_idx = 1;
  bool in_domain = false;
  double score = 0;
  for (size_t i = 0; i < classes.size(); ++i) {
    if (reset_points[reset_idx] == i) {
      if (in_domain) {
        scores.push_back(score);
        score = 0;
        in_domain = false;
      }
      ++reset_idx;
    }
    if (classes[i]) {
      for (size_t r = 0; r < n_replicates ; ++r)
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
build_domains(const vector<SimpleGenomicRegion> &bins,
              const vector<size_t> &reset_points,
              const vector<bool> &classes,
              vector<GenomicRegion> &domains) {

  size_t n_bins = 0, reset_idx = 1, prev_end = 0;
  bool in_domain = false;
  for (size_t i = 0; i < classes.size(); ++i) {
    if (reset_points[reset_idx] == i) {
      if (in_domain) {
        domains.back().set_end(prev_end);
        domains.back().set_score(n_bins);
        n_bins = 0;
        in_domain = false;
      }
      ++reset_idx;
    }
    if (classes[i]) {
      if (!in_domain) {
        domains.push_back(GenomicRegion(bins[i]));
        in_domain = true;
      }
      ++n_bins;
    }
    else if (in_domain) {
      domains.back().set_end(prev_end);
      domains.back().set_score(n_bins);
      n_bins = 0;
      in_domain = false;
    }
    prev_end = bins[i].get_end();
  }
}


//Modified to take multiple replicates
template <class T, class U> static void
separate_regions(const size_t desert_size,
                 vector<vector<SimpleGenomicRegion> > &bins,
                 vector<vector<T> > &meth, vector<vector<U>> &reads,
                 vector<size_t> &reset_points,
                 vector<size_t> &dists_btwn_bins) {
  const size_t n_replicates = bins.size();

  // eliminate the zero-read cpg sites if no coverage in any rep
  size_t end_coord_of_prev = 0;
  size_t j = 0;
  for (size_t i = 0; i < bins[0].size(); ++i) {
    bool all_empty = true;
    size_t rep_idx = 0;
    while (all_empty && rep_idx < n_replicates) {
      if (reads[rep_idx][i] == 0) ++rep_idx;
      else all_empty = false;
    }
    if (!all_empty) {
      dists_btwn_bins.push_back(bins[0][i].get_start() - end_coord_of_prev);
      end_coord_of_prev = bins[0][i].get_end();
      for (size_t r = 0; r < n_replicates; ++r) {
        bins[r][j] = bins[r][i];
        meth[r][j] = meth[r][i];
        reads[r][j] = reads[r][i];
      }
      ++j;
    }
  }

  for (size_t r = 0; r < n_replicates; ++r) {
    bins[r].resize(j);
    meth[r].resize(j);
    reads[r].resize(j);
  }

  // segregate bins
  size_t prev_cpg = 0;
  for (size_t i = 0; i < bins[0].size(); ++i) {
    const size_t dist = (i > 0 && bins[0][i].same_chrom(bins[0][i - 1])) ?
      bins[0][i].get_start() - prev_cpg : num_lim<size_t>::max();
    if (dist > desert_size)
      reset_points.push_back(i);
    prev_cpg = bins[0][i].get_start();
  }
  assert(std::size(reset_points) > 0);
  reset_points.push_back(bins[0].size());
}


static void
shuffle_bins(const size_t rng_seed,
             const TwoStateHMM &hmm,
             vector<vector<pair<double, double> > > meth,
             const vector<size_t> &reset_points,
             const vector<double> &start_trans,
             const vector<vector<double> > &trans,
             const vector<double> &end_trans,
             const vector<double> &fg_alpha, const vector<double> &fg_beta,
             const vector<double> &bg_alpha, const vector<double> &bg_beta,
             vector<double> &domain_scores,
             vector<bool> &array_status) {

  size_t n_replicates = meth.size();

  auto eng = std::default_random_engine(rng_seed);
  for (size_t r =0 ; r < n_replicates; ++r)
    std::shuffle(begin(meth[r]), end(meth[r]), eng);

  vector<bool> classes;
  vector<double> scores;
  hmm.PosteriorDecoding_rep(meth, reset_points, start_trans, trans,
                            end_trans, fg_alpha, fg_beta, bg_alpha,
                            bg_beta, classes, scores, array_status);
  get_domain_scores(classes, meth, reset_points, domain_scores);
  sort(begin(domain_scores), end(domain_scores));
}


static void
assign_p_values(const vector<double> &random_scores,
                const vector<double> &observed_scores,
                vector<double> &p_values) {
  const double n_randoms = random_scores.size() == 0 ? 1 : random_scores.size();
  for (auto scr : observed_scores) {
    const auto scr_itr =
      upper_bound(cbegin(random_scores), cend(random_scores), scr);
    p_values.push_back(std::distance(scr_itr, cend(random_scores)) / n_randoms);
  }
}


static void
read_params_file(const bool verbose,
                 const string &params_file,
                 double &fg_alpha,
                 double &fg_beta,
                 double &bg_alpha,
                 double &bg_beta,
                 vector<double> &start_trans,
                 vector<vector<double> > &trans,
                 vector<double> &end_trans,
                 double &fdr_cutoff) {
  string jnk;
  std::ifstream in(params_file.c_str());

  in >> jnk >> fg_alpha
     >> jnk >> fg_beta
     >> jnk >> bg_alpha
     >> jnk >> bg_beta
     >> jnk >> start_trans[0]
     >> jnk >> start_trans[1]
     >> jnk >> trans[0][0]
     >> jnk >> trans[0][1]
     >> jnk >> trans[1][0]
     >> jnk >> trans[1][1]
     >> jnk >> end_trans[0]
     >> jnk >> end_trans[1]
     >> jnk >> fdr_cutoff;

  if (verbose)
    cerr << "Read in params from " << params_file << endl
         << "FG_ALPHA\t" << fg_alpha << endl
         << "FG_BETA\t" << fg_beta << endl
         << "BG_ALPHA\t" << bg_alpha << endl
         << "BG_BETA\t" << bg_beta << endl
         << "S_F\t" << start_trans[0] << endl
         << "S_B\t" << start_trans[1] << endl
         << "F_F\t" << trans[0][0] << endl
         << "F_B\t" << trans[0][1] << endl
         << "B_F\t" << trans[1][0] << endl
         << "B_B\t" << trans[1][1] << endl
         << "F_E\t" << end_trans[0] << endl
         << "B_E\t" << end_trans[1] << endl
         << "FDR_CUTOFF\t" << fdr_cutoff << endl;
}


static void
write_posteriors_file(const string &posteriors_file,
                      const vector<vector<SimpleGenomicRegion> > &bins,
                      const vector<double> &scores) {
  static const size_t decimal_precision = 10;

  ofstream out(posteriors_file);
  out.precision(decimal_precision);
  for (size_t r = 0; r < scores.size(); ++r)
    out << bins[0][r] << '\t' << scores[r] << endl;
}


static void
write_params_file(const string &outfile,
                  const vector<double> &fg_alpha,
                  const vector<double> &fg_beta,
                  const vector<double> &bg_alpha,
                  const vector<double> &bg_beta,
                  const vector<double> &start_trans,
                  const vector<vector<double> > &trans,
                  const vector<double> &end_trans) {
  static const size_t decimal_precision = 30;
  ofstream out(outfile);
  out.precision(decimal_precision);
  for (size_t r =0; r < fg_alpha.size(); ++r)
    out << "FG_ALPHA_" << r+1 << "\t" << std::setw(14) << fg_alpha[r] << "\t"
        << "FG_BETA_" << r+1 << "\t" << std::setw(14) << fg_beta[r] << "\t"
        << "BG_ALPHA_" << r+1 << "\t" << std::setw(14) << bg_alpha[r] << "\t"
        << "BG_BETA_" << r+1 << "\t" << std::setw(14) << bg_beta[r] << endl;

  out << "S_F\t" << start_trans[0] << endl
      << "S_B\t" << start_trans[1] << endl
      << "F_F\t" << trans[0][0] << endl
      << "F_B\t" << trans[0][1] << endl
      << "B_F\t" << trans[1][0] << endl
      << "B_B\t" << trans[1][1] << endl
      << "F_E\t" << end_trans[0] << endl
      << "B_E\t" << end_trans[1] << endl;
}


static bool
check_if_array_data(const string &infile) {

  bgzf_file in(infile, "r");
  if (!in) throw std::runtime_error("bad file: " + infile);

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
load_array_data(const size_t bin_size,
                const string &cpgs_file,
                vector<SimpleGenomicRegion> &bins,
                vector<pair<double, double> > &meth,
                vector<size_t> &reads) {
  // MAGIC. GS: minimum value for array?
  static const double meth_min = 1.0e-2;

  bgzf_file in(cpgs_file, "r");
  if (!in) throw std::runtime_error("bad sites file: " + cpgs_file);

  if (get_has_counts_header(cpgs_file))
    skip_counts_header(in);

  string curr_chrom;
  size_t prev_pos = 0ul, curr_pos = 0ul;
  double array_meth_bin = 0.0;
  double num_probes_in_bin = 0.0;

  MSite site;
  while (read_site(in, site)) {
    // TODO: MN: I think that the block below should be placed later
    // in this scope. At this location, the methylation level of the
    // first site in a new chrom is contributed to the last bin of the
    // previous chrom.
    if (site.n_reads > 0) { // its covered by a probe
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
          throw runtime_error("CpGs not sorted in file \""
                              + cpgs_file + "\"");
        bins.push_back(SimpleGenomicRegion(curr_chrom, curr_pos,
                                           prev_pos + 1));
        meth.push_back(make_pair(array_meth_bin, num_probes_in_bin));
        if (num_probes_in_bin > 0) reads.push_back(1);
        else reads.push_back(0);
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
      bins.push_back(SimpleGenomicRegion(curr_chrom, curr_pos,
                                         curr_pos + bin_size));
      meth.push_back(make_pair(array_meth_bin, num_probes_in_bin));
      (num_probes_in_bin > 0) ? reads.push_back(1) : reads.push_back(0);

      array_meth_bin = 0.0;
      num_probes_in_bin = 0.0;
      curr_pos += bin_size;
      while (curr_pos + bin_size < site.pos) {
        bins.push_back(SimpleGenomicRegion(curr_chrom, curr_pos,
                                           curr_pos + bin_size));
        reads.push_back(0);
        meth.push_back(make_pair(0.0, 0.0));
        curr_pos += bin_size;
      }
    }
  }

  if (site.meth != -1 ) { // its covered by a probe
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
    bins.push_back(SimpleGenomicRegion(curr_chrom, curr_pos,
                                       prev_pos + 1));
    meth.push_back(make_pair(array_meth_bin, num_probes_in_bin));
    if (num_probes_in_bin > 0)
      reads.push_back(1);
    else reads.push_back(0);
  }
}


static void
load_wgbs_data(const size_t bin_size, const string &cpgs_file,
               vector<SimpleGenomicRegion> &bins,
               vector<pair<double, double> > &meth,
               vector<size_t> &reads) {
  reads.clear(); // for safety
  meth.clear();
  bins.clear();

  // ADS: loading data each iteration should be put outside loop
  bgzf_file in(cpgs_file, "r");
  if (!in) throw runtime_error("bad sites file: " + cpgs_file);

  if (get_has_counts_header(cpgs_file))
    skip_counts_header(in);

  // keep track of the chroms we've seen
  string curr_chrom;
  std::unordered_set<string> chroms_seen;

  MSite site;
  size_t prev_pos = 0ul;
  size_t sites_in_bin = 0ul;

  while (read_site(in, site)) {
    if (curr_chrom != site.chrom) { // handle change of chrom
      if (sites_in_bin > 0) bins.back().set_end(prev_pos);
      if (chroms_seen.find(site.chrom) != end(chroms_seen))
        throw runtime_error("sites not sorted");
      chroms_seen.insert(site.chrom);
      curr_chrom = site.chrom;
      reads.push_back(0);
      meth.push_back(make_pair(0.0, 0.0));
      bins.push_back(SimpleGenomicRegion(site.chrom, 0, bin_size));
      sites_in_bin = 0;
    }
    prev_pos = site.pos;
    if (site.pos < bins.back().get_start())
      throw runtime_error("sites not sorted");
    while (bins.back().get_end() < site.pos) {
      sites_in_bin = 0;
      reads.push_back(0);
      meth.push_back(make_pair(0.0, 0.0));
      bins.push_back(SimpleGenomicRegion(site.chrom, bins.back().get_end(),
                                         bins.back().get_end() + bin_size));
    }
    reads.back() += site.n_reads;
    meth.back().first += site.n_meth();
    meth.back().second += site.n_unmeth();
    sites_in_bin++;
  }
  if (sites_in_bin > 0) bins.back().set_end(prev_pos);
}


static void
remove_empty_bins_at_chrom_start(vector<SimpleGenomicRegion> &bins,
                                 vector<pair<double, double> > &meth,
                                 vector<size_t> &reads) {
  bool chrom_start = true;
  size_t j = 0;
  string prev_chrom = "";
  for (size_t i = 0; i < bins.size(); ++i) {
    if (bins[i].get_chrom() != prev_chrom) {
      chrom_start = true;
      prev_chrom = bins[i].get_chrom();
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
  bins.erase(begin(bins) + j, end(bins));
  meth.erase(begin(meth) + j, end(meth));
  reads.erase(begin(reads) + j, end(reads));
}


static void
load_read_counts(const string &cpgs_file, const size_t bin_size,
                 vector<size_t> &reads) {
  reads.clear(); // for safety

  // ADS: loading data each iteration should be put outside loop
  bgzf_file in(cpgs_file, "r");
  if (!in) throw runtime_error("bad methcounts file: " + cpgs_file);

  if (get_has_counts_header(cpgs_file))
    skip_counts_header(in);

  // keep track of where we are and what we've seen
  size_t bin_start = 0ul;
  string curr_chrom;
  std::unordered_set<string> chroms_seen;

  MSite site;
  while (read_site(in, site)) {
    if (curr_chrom != site.chrom) { // handle change of chrom
      if (chroms_seen.find(site.chrom) != end(chroms_seen))
        throw runtime_error("sites not sorted");
      chroms_seen.insert(site.chrom);
      bin_start = 0;
      curr_chrom = site.chrom;
      reads.push_back(0);
    }
    if (site.pos < bin_start)
      throw runtime_error("sites not sorted");
    for (; bin_start + bin_size < site.pos; bin_start += bin_size)
      reads.push_back(0);
    reads.back() += site.n_reads;
  }
}


static double
good_bins_frac(const vector<size_t> &cumulative, const size_t min_bin_size,
               const size_t bin_size, const size_t min_cov_to_pass) {

  // make sure the target bin size is a multiple of the minimum so we
  // have the resolution to construct the new bins
  assert(bin_size % min_bin_size == 0);

  // the step size corresponds to the number of minium sized bins that
  // would make up a new bin of the target size
  const size_t step_size = bin_size/min_bin_size;

  size_t passing_bins = 0, covered_bins = 0;

  size_t prev_total = 0;
  for (size_t i = 0; i + step_size < cumulative.size(); i += step_size) {
    const size_t curr_cumulative = cumulative[i + step_size];
    const size_t bin_count = curr_cumulative - prev_total;
    covered_bins += (bin_count > 0);
    passing_bins += (bin_count >= min_cov_to_pass);
    prev_total = curr_cumulative;
  }

  // check if there is a leftover bin at the end
  if (cumulative.size() % step_size != 0) {
    const size_t bin_count = cumulative.back() - prev_total;
    covered_bins += (bin_count > 0);
    passing_bins += (bin_count >= min_cov_to_pass);
  }

  return static_cast<double>(passing_bins)/std::max(1ul, covered_bins);
}

static size_t
get_min_reads_for_confidence(const double conf_level) {
  // ADS: value of 0.5 below important; this is where the CI is widest
  static const double fixed_phat = 0.5;
  size_t n_reads = 0;
  double lower = 0.0, upper = 1.0;
  // ADS: should be doubling first, followed by bisection
  while (1.0 - conf_level < upper - lower) {
    ++n_reads;
    wilson_ci_for_binomial(1.0 - conf_level, n_reads,
                           fixed_phat, lower, upper);
  }
  return n_reads;
}


// ADS: this function will return num_lim<size_t>::max() if the
// fraction of "good" bins is zero for all attempted bin sizes.
static size_t
binsize_selection(const size_t resolution, const size_t min_bin_sz,
                  const size_t max_bin_sz, const double conf_level,
                  const double min_frac_passed, const string &cpgs_file) {

  const size_t min_cov_to_pass = get_min_reads_for_confidence(conf_level);

  vector<size_t> reads;
  load_read_counts(cpgs_file, resolution, reads);

  std::partial_sum(begin(reads), end(reads), begin(reads));

  double frac_passed = 0.0;
  size_t bin_size = min_bin_sz;

  while (bin_size < max_bin_sz && frac_passed < min_frac_passed) {
    frac_passed = good_bins_frac(reads, resolution, bin_size, min_cov_to_pass);
    if (frac_passed < min_frac_passed)
      bin_size += resolution;
  }
  return frac_passed < min_frac_passed ? num_lim<size_t>::max() : bin_size;
}


static void
load_bins(const size_t bin_size,
          const string &cpgs_file,
          vector<SimpleGenomicRegion> &bins,
          vector<pair<double, double> > &meth,
          vector<size_t> &reads, vector<bool> &array_status) {

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
get_union_of_bins(const vector<vector<SimpleGenomicRegion> > &orig,
                  vector<SimpleGenomicRegion> &bins) {

  // flatten the set of sorted bins
  bins.clear();
  for (auto &&i: orig)
    bins.insert(end(bins), begin(i), end(i));

  // merge each sorted interval of bins
  const auto first = begin(bins);
  auto middle = begin(bins);
  for (size_t i = 1; i < orig.size(); ++i) {
    middle += orig[i-1].size();
    std::inplace_merge(first, middle, middle + orig[i].size());
  }
  // ensure unique bins
  bins.erase(unique(begin(bins), end(bins)), end(bins));
  bins.shrink_to_fit();

  // make sure all bins are aligned at same boundaries
  for (size_t i = 1; i < bins.size(); ++i)
    if (bins[i-1].overlaps(bins[i]))
      throw std::runtime_error("bins from reps not aligned");
}


static void
add_missing_bins(const vector<SimpleGenomicRegion> &all_bins,
                 vector<SimpleGenomicRegion> &bins,
                 vector<pair<double, double>> &meth) {

  const size_t n_bins = all_bins.size();
  vector<pair<double, double>> tmp_meth(n_bins);

  size_t j = 0; // assume j range no larger than i range
  for (size_t i = 0; i < n_bins; ++i) {
    if (all_bins[i] == bins[j])
      tmp_meth[i] = meth[j++];
    else
      tmp_meth[i] = {0.0, 0.0};
  }
  std::swap(meth, tmp_meth);
  bins = all_bins;
}


static void
write_empty_summary(const string &summary_file) {
  if (!summary_file.empty()) {
    ofstream summary_out(summary_file);
    if (!summary_out)
      throw runtime_error("failed to open: " + summary_file);
    summary_out << pmd_summary({}).tostring() << endl;
  }
}


int
main_pmd(int argc, char *argv[]) {
  try {

    static const size_t min_observations_for_inference = 100;
    static const size_t max_bin_size = 500000;
    static const size_t min_bin_size = 1000;
    size_t resolution = 500;

    const char* sep = ",";
    string outfile;

    size_t rng_seed = 408;

    bool DEBUG = false;
    size_t desert_size = 5000;
    size_t bin_size = 1000;
    size_t max_iterations = 10;
    // run mode flags
    bool verbose = false;
    bool ARRAY_MODE = false;
    bool fixed_bin_size = false;

    // MAGIC: corrections for small values (not parameters):
    static const double tolerance = 1e-5;
    static const double min_prob  = 1e-10;

    string summary_file;

    string params_in_files;
    string params_out_file;
    string posteriors_out_prefix;


    const string description =
      "Identify PMDs in methylomes. Methylation must be provided in the \
      methcounts file format (chrom, position, strand, context, \
      methylation, reads). See the methcounts documentation for \
      details. This program assumes only data at CpG sites and that \
      strands are collapsed so only the positive site appears in the \
      file, but reads counts are from both strands.";

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), description,
                           "<methcount-files>");
    opt_parse.add_opt("out", 'o', "output file (default: stdout)",
                      false, outfile);
    opt_parse.add_opt("desert", 'd', "max dist between bins with data in PMD",
                      false, desert_size);
    opt_parse.add_opt("fixedbin", 'f', "Fixed bin size", false, fixed_bin_size);
    opt_parse.add_opt("bin", 'b', "Starting bin size", false, bin_size);
    opt_parse.add_opt("arraymode",'a', "All samples are array",
                      false, ARRAY_MODE);
    opt_parse.add_opt("itr", 'i', "max iterations", false, max_iterations);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, verbose);
    opt_parse.add_opt("debug", 'D', "print more run info", false, DEBUG);
    opt_parse.add_opt("params-in", 'P', "HMM parameter files for "
                      "individual methylomes (separated with comma)",
                      false, params_in_files);
    opt_parse.add_opt("posteriors-out", 'r',
                      "write out posterior probabilities in methcounts format",
                      false, posteriors_out_prefix);
    opt_parse.add_opt("summary", 'S',
                      "write summary output here",
                      false, summary_file);
    opt_parse.add_opt("params-out", 'p', "write HMM parameters to this file",
                      false, params_out_file);
    opt_parse.add_opt("seed", 's', "specify random seed",
                      false, rng_seed);

    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl
           << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    if (leftover_args.empty()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    /****************** END COMMAND LINE OPTIONS *****************/

    resolution = min(bin_size, resolution);

    vector<string> cpgs_file = leftover_args;
    vector<string> params_in_file;
    if (!params_in_files.empty()) {
      params_in_file = smithlab::split(params_in_files, sep, false);
      assert(cpgs_file.size() == params_in_file.size());
    }

    const size_t n_replicates = cpgs_file.size();
    for (auto &filename : cpgs_file)
      if (!is_msite_file(filename))
        throw runtime_error("malformed counts file: " + filename);

    bool insufficient_data = false; // ADS: this is used now to detect
                                    // when the counts files have
                                    // lines for CpG sites, but no
                                    // counts.

    // Sanity checks input file format and dynamically selects bin
    // size from WGBS samples.
    if (!fixed_bin_size && !ARRAY_MODE) {
      if (verbose)
        cerr << "[DYNAMICALLY SELECTING BIN SIZE]" << endl;
      double confidence_interval = 0.80;
      double prop_accept = 0.80;
      for (size_t i = 0; i < n_replicates && !insufficient_data; ++i) {
        const bool arrayData = check_if_array_data(cpgs_file[i]);
        if (!arrayData) {
          bin_size = binsize_selection(resolution, min_bin_size, max_bin_size,
                                       confidence_interval, prop_accept,
                                       cpgs_file[i]);
          if (bin_size == num_lim<size_t>::max())
            insufficient_data = true;
          desert_size = 5*bin_size; // TODO: explore extrapolation number
        }
        else {
          // same as the parameters below
          bin_size = 1000;
          desert_size = 200000;
        }
      }
    }
    else if (ARRAY_MODE) {
      bin_size = 1000;    // MAGIC NUMBERS FROM PAPER
      desert_size = 200000; // PERFORM WITH HIGHEST JACCARD INDEX TO WGBS
    }
    else {
      desert_size = max(desert_size, bin_size);
    }

    if (insufficient_data) {
      // ADS: first check for insufficient data; another is needed if
      // fixed bin size is used
      if (verbose) cerr << "EXITING: INSUFFICIENT DATA" << endl;
      if (!summary_file.empty()) write_empty_summary(summary_file);
      return EXIT_SUCCESS;
    }

    if (verbose)
      cerr << "[READING IN AT BIN SIZE " << bin_size << "]" << endl;

    // separate the regions by chrom and by desert
    vector<vector<SimpleGenomicRegion> > bins(n_replicates);
    vector<vector<pair<double, double> > > meth(n_replicates);
    vector<vector<size_t> > reads(n_replicates);
    vector<bool> array_status;

    for (size_t i = 0; i < n_replicates && !insufficient_data; ++i) {
      if (verbose)
        cerr << "[READING CPGS AND METH PROPS] from " << cpgs_file[i] << endl;

      load_bins(bin_size, cpgs_file[i], bins[i], meth[i],
                reads[i], array_status);
      const double total_observations =
        accumulate(begin(reads[i]), end(reads[i]), 0);
      if (total_observations <= num_lim<double>::min())
        insufficient_data = true;
      if (verbose)
        cerr << "TOTAL BINS: " << bins[i].size() << endl
             << "MEAN COVERAGE: "
             << total_observations/std::max(1ul, reads[i].size())
             << endl;
    }

    if (insufficient_data) {
      // ADS: second check for insufficient data; another is needed if
      // filtered number of bins is too few
      if (verbose) cerr << "EXITING: INSUFFICIENT DATA" << endl;
      if (!summary_file.empty()) write_empty_summary(summary_file);
      return EXIT_SUCCESS;
    }

    if (n_replicates > 1) {
      bool need_to_adjust_bins = false;
      const size_t n_bins_0 = bins[0].size();
      for (size_t i = 1; i < n_replicates && !need_to_adjust_bins; ++i)
        if (bins[i].size() != n_bins_0)
          need_to_adjust_bins = true;

      if (need_to_adjust_bins) {
        vector<SimpleGenomicRegion> all_bins;
        get_union_of_bins(bins, all_bins);
        for (size_t i = 0; i < bins.size(); ++i)
          add_missing_bins(all_bins, bins[i], meth[i]);
      }
    }

    // separate regions by chrom and desert; eliminate isolated Bins
    vector<size_t> reset_points;
    vector<size_t> dists_btwn_bins;
    if (verbose)
      cerr << "[separating by CpG desert]" << endl;
    separate_regions(desert_size, bins, meth, reads,
                     reset_points, dists_btwn_bins);
    if (size(bins[0]) < min_observations_for_inference)
      insufficient_data = true;

    if (insufficient_data) {
      // ADS: final check for sufficient data failed; too few bins
      // after filtering
      if (verbose) cerr << "EXITING: INSUFFICIENT DATA" << endl;
      if (!summary_file.empty()) write_empty_summary(summary_file);
      return EXIT_SUCCESS;
    }

    if (verbose)
      cerr << "bins retained: " << std::size(bins[0]) << endl
           << "number of distances between: " << std::size(dists_btwn_bins) << endl
           << "deserts removed: " << size(reset_points) - 2 << endl;

    /****************** Read in params *****************/
    vector<double> start_trans(2, 0.5), end_trans(2, 1e-10);
    vector<vector<double> > trans(2, vector<double>(2, 0.01));
    trans[0][0] = trans[1][1] = 0.99;
    const TwoStateHMM hmm(min_prob, tolerance, max_iterations, verbose, DEBUG);
    vector<double> reps_fg_alpha(n_replicates, 0.05);
    vector<double> reps_fg_beta(n_replicates, 0.95);
    vector<double> reps_bg_alpha(n_replicates, 0.95);
    vector<double> reps_bg_beta(n_replicates, 0.05);
    double score_cutoff_for_fdr = num_lim<double>::max();

    if (!params_in_file.empty()) {
      // read parameters files
      for (size_t i= 0; i < n_replicates; ++i)
        read_params_file(verbose, params_in_file[i], reps_fg_alpha[i],
                         reps_fg_beta[i], reps_bg_alpha[i], reps_bg_beta[i],
                         start_trans, trans, end_trans, score_cutoff_for_fdr);
    }

    // train model (default behavior; not done when params supplied)
    if (max_iterations > 0)
      hmm.BaumWelchTraining_rep(meth, reset_points,
                                start_trans, trans, end_trans,
                                reps_fg_alpha, reps_fg_beta,
                                reps_bg_alpha, reps_bg_beta, array_status);

    if (!params_out_file.empty()) {
      // write all the HMM parameters
      write_params_file(params_out_file,
                        reps_fg_alpha, reps_fg_beta,
                        reps_bg_alpha, reps_bg_beta,
                        start_trans, trans, end_trans);
    }

    /***********************************/

    if (!posteriors_out_prefix.empty()) {
      vector<double> into_scores;
      hmm.TransitionPosteriors_rep(meth, reset_points, start_trans, trans,
                                   end_trans, reps_fg_alpha, reps_fg_beta,
                                   reps_bg_alpha, reps_bg_beta, array_status,
                                   2, into_scores);
      write_posteriors_file(posteriors_out_prefix + ".intoTrans", bins,
                            into_scores);
      vector<double> outof_scores;
      hmm.TransitionPosteriors_rep(meth, reset_points, start_trans, trans,
                                   end_trans, reps_fg_alpha, reps_fg_beta,
                                   reps_bg_alpha, reps_bg_beta, array_status,
                                   1, outof_scores);
      write_posteriors_file(posteriors_out_prefix + ".outofTrans", bins,
                            outof_scores);
    }

    vector<bool> classes;
    vector<double> scores;
    hmm.PosteriorDecoding_rep(meth, reset_points, start_trans, trans, end_trans,
                              reps_fg_alpha, reps_fg_beta,
                              reps_bg_alpha, reps_bg_beta, classes,
                              scores, array_status);

    if (!posteriors_out_prefix.empty())
      write_posteriors_file(posteriors_out_prefix + ".posteriors", bins, scores);

    vector<double> domain_scores;
    get_domain_scores(classes, meth, reset_points, domain_scores);

    if (verbose)
      cerr << "[RANDOMIZING SCORES FOR FDR]" << endl;

    vector<double> random_scores;
    shuffle_bins(rng_seed, hmm, meth, reset_points, start_trans, trans,
                 end_trans, reps_fg_alpha, reps_fg_beta, reps_bg_alpha,
                 reps_bg_beta, random_scores, array_status);

    vector<double> p_values;
    assign_p_values(random_scores, domain_scores, p_values);

    if (score_cutoff_for_fdr == num_lim<double>::max() &&
        !p_values.empty())
      score_cutoff_for_fdr = get_score_cutoff_for_fdr(p_values, 0.01);

    if (!params_out_file.empty()) {
      ofstream out(params_out_file, std::ios::app);
      out << "FDR_CUTOFF\t"
          << std::setprecision(30) << score_cutoff_for_fdr << endl;
    }

    vector<GenomicRegion> domains;
    build_domains(bins[0], reset_points, classes, domains);

    size_t good_pmd_count = 0;
    vector<GenomicRegion> good_domains;
    for (size_t i = 0; i < domains.size(); ++i)
      if (p_values[i] < score_cutoff_for_fdr) {
        good_domains.push_back(domains[i]);
        good_domains.back().set_name("PMD" + to_string(good_pmd_count++));
      }

    optimize_boundaries(bin_size, cpgs_file, good_domains, array_status,
                        reps_fg_alpha, reps_fg_beta,
                        reps_bg_alpha, reps_bg_beta);

    if (!summary_file.empty()) {
      ofstream summary_out(summary_file);
      if (!summary_out) throw runtime_error("failed to open: " + summary_file);
      summary_out << pmd_summary(good_domains).tostring() << endl;
    }

    ofstream of;
    if (!outfile.empty()) of.open(outfile);
    ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

    copy(begin(good_domains), end(good_domains),
         std::ostream_iterator<GenomicRegion>(out, "\n"));
  }
  catch (const runtime_error &e) {
    cerr << "ERROR:\t" << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
