/*
  Copyright (C) 2011-2023 University of Southern California
                          Andrew D. Smith and Song Qiang

  Authors: Andrew D. Smith, Song Qiang

  This file is part of dnmtools.

  dnmtools is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  dnmtools is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#ifndef THREE_STATE_HMM_HPP
#define THREE_STATE_HMM_HPP

#include <memory>
#include <utility>
#include <string>
#include <vector>

#include "smithlab_utils.hpp"
#include "BetaBin.hpp"

enum STATE_LABELS {hypo, HYPER, HYPO};

struct Triplet {double hypo, HYPER, HYPO;};

class ThreeStateHMM {
public:

  ThreeStateHMM(std::vector<std::pair<double, double>> &obs,
                const std::vector<size_t> &res,
                const double tol, const size_t max_itr, const bool v);

  void
  set_parameters(const betabin & hypo_em,
                 const betabin & HYPER_em,
                 const betabin & HYPO_em,
                 const std::vector<std::vector<double>> &tr);

  void
  get_parameters(betabin & hypo_em,
                 betabin & HYPER_em,
                 betabin & HYPO_em,
                 std::vector<std::vector<double> > &tr) const;

  double BaumWelchTraining();

  double PosteriorDecoding();

  double ViterbiDecoding();

  void get_state_posteriors(std::vector<Triplet> &scores) const;

  void get_classes(std::vector<STATE_LABELS>  &classes) const;

  // private:

  //////////// methods ////////////
  double single_iteration();
  double forward_algorithm(const size_t start, const size_t end);
  double backward_algorithm(const size_t start, const size_t end);
  double hypo_segment_log_likelihood(const size_t start, const size_t end);
  double HYPER_segment_log_likelihood(const size_t start, const size_t end);
  double HYPO_segment_log_likelihood(const size_t start, const size_t end);

  void estimate_state_posterior(const size_t start, const size_t end);
  void estimate_posterior_trans_prob(const size_t start, const size_t end);
  void estimate_parameters();
  void update_observation_likelihood();

  double ViterbiDecoding(const size_t start, const size_t end);

  ////////   data   ////////
  std::vector<std::pair<double, double>> observations;
  std::vector<size_t> reset_points;
  std::vector<double> meth_lp, unmeth_lp;
  std::vector<double> hypo_log_likelihood, HYPER_log_likelihood, HYPO_log_likelihood;

  //  HMM internal data
  betabin hypo_emission, HYPER_emission, HYPO_emission;

  Triplet lp_start, lp_end;
  std::vector<std::vector<double> > trans;

  std::vector<Triplet> forward;
  std::vector<Triplet> backward;
  std::vector<double> hypo_posteriors, HYPER_posteriors, HYPO_posteriors;
  std::vector<double> hypo_hypo, hypo_HYPER,
    HYPER_hypo, HYPER_HYPER, HYPER_HYPO,
    HYPO_HYPER, HYPO_HYPO;

  // result
  std::vector<STATE_LABELS> classes;
  std::vector<Triplet> state_posteriors;

  // parameters
  double tolerance;
  size_t max_iterations;
  bool VERBOSE;
};

#endif
