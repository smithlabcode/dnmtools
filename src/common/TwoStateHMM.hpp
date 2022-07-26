/*
  Copyright (C) 2019-2022 Andrew D. Smith
  Author: Andrew D. Smith

  This is free software; you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This software is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.
*/

#ifndef TWO_STATE_HMM_HPP
#define TWO_STATE_HMM_HPP

#include <memory>
#include <vector>

struct TwoStateBetaBin;

class TwoStateHMM {
public:

  TwoStateHMM(const double tol, const size_t max_itr, const bool v) :
    tolerance(tol), max_iterations(max_itr), VERBOSE(v) {}

  double
  ViterbiDecoding(const std::vector<std::pair<double, double> > &values,
                  const std::vector<size_t> &reset_points,
                  const double f_to_b_trans, const double b_to_f_trans,
                  const double fg_alpha, const double fg_beta,
                  const double bg_alpha, const double bg_beta,
                  std::vector<bool> &ml_classes) const;


  double
  BaumWelchTraining(const std::vector<std::pair<double, double> > &values,
                    const std::vector<size_t> &reset_points,
                    double &f_to_b_trans, double &b_to_f_trans,
                    double &fg_alpha, double &fg_beta,
                    double &bg_alpha, double &bg_beta) const;

  double
  PosteriorDecoding(const std::vector<std::pair<double, double> > &values,
                    const std::vector<size_t> &reset_points,
                    const double f_to_b_trans, const double b_to_f_trans,
                    const double fg_alpha, const double fg_beta,
                    const double bg_alpha, const double bg_beta,
                    std::vector<bool> &classes,
                    std::vector<double> &llr_scores) const;

  void
  PosteriorScores(const std::vector<std::pair<double, double> > &values,
                  const std::vector<size_t> &reset_points,
                  const double f_to_b_trans, const double b_to_f_trans,
                  const double fg_alpha, const double fg_beta,
                  const double bg_alpha, const double bg_beta,
                  const bool class_id,
                  std::vector<double> &llr_scores) const;

  void
  TransitionPosteriors(const std::vector<std::pair<double, double> > &values,
                       const std::vector<size_t> &reset_points,
                       const double f_to_b_trans, const double b_to_f_trans,
                       const double fg_alpha, const double fg_beta,
                       const double bg_alpha, const double bg_beta,
                       const size_t transition,
                       std::vector<double> &scores) const;

  // FOR MULTIPLE REPLICATES
  double
  BaumWelchTraining(const std::vector<std::vector<std::pair<double, double> > > &values,
                    const std::vector<size_t> &reset_points,
                    double &f_to_b_trans, double &b_to_f_trans,
                    std::vector<double> &fg_alpha,
                    std::vector<double> &fg_beta,
                    std::vector<double> &bg_alpha,
                    std::vector<double> &bg_beta) const;

  double
  PosteriorDecoding(const std::vector<std::vector<std::pair<double, double> > > &values,
                    const std::vector<size_t> &reset_points,
                    const double f_to_b_trans, const double b_to_f_trans,
                    const std::vector<double> &fg_alpha,
                    const std::vector<double> &fg_beta,
                    const std::vector<double> &bg_alpha,
                    const std::vector<double> &bg_beta,
                    std::vector<bool> &classes,
                    std::vector<double> &llr_scores) const;

  void
  PosteriorScores(const std::vector<std::vector<std::pair<double, double> > > &values,
                  const std::vector<size_t> &reset_points,
                  const double f_to_b_trans, const double b_to_f_trans,
                  const std::vector<double> &fg_alpha,
                  const std::vector<double> &fg_beta,
                  const std::vector<double> &bg_alpha,
                  const std::vector<double> &bg_beta,
                  const bool &fg_class,
                  std::vector<double> &llr_scores) const;


private:

  double
  ViterbiDecoding(const std::vector<std::pair<double, double> > &values,
                  const std::vector<size_t> &reset_points,
                  const double p_fb, const double p_bf,
                  const TwoStateBetaBin &fg_distro, const TwoStateBetaBin &bg_distro,
                  std::vector<bool> &ml_classes) const;

  double
  BaumWelchTraining(const std::vector<std::pair<double, double> > &values,
                    const std::vector<size_t> &reset_points,
                    double &p_fb, double &p_bf,
                    TwoStateBetaBin &fg_distro, TwoStateBetaBin &bg_distro) const;

  double
  PosteriorDecoding(const std::vector<std::pair<double, double> > &values,
                    const std::vector<size_t> &reset_points,
                    const double p_fb, const double p_bf,
                    const TwoStateBetaBin &fg_distro,
                    const TwoStateBetaBin &bg_distro,
                    std::vector<bool> &classes,
                    std::vector<double> &llr_scores) const;

  void
  PosteriorScores(const std::vector<std::pair<double, double> > &values,
                  const std::vector<size_t> &reset_points,
                  const double p_fb, const double p_bf,
                  const TwoStateBetaBin &fg_distro, const TwoStateBetaBin &bg_distro,
                  const bool class_id,
                  std::vector<double> &llr_scores) const;

  void
  TransitionPosteriors(const std::vector<std::pair<double, double> > &values,
                       const std::vector<size_t> &reset_points,
                       const double p_fb, const double p_bf,
                       const TwoStateBetaBin &fg_distro, const TwoStateBetaBin &bg_distro,
                       const size_t transition,
                       std::vector<double> &scores) const;

  // FOR MULTIPLE REPLICATES

  double
  BaumWelchTraining(const std::vector<std::vector<std::pair<double, double> > > &values,
                    const std::vector<size_t> &reset_points,
                    double &p_fb, double &p_bf,
                    std::vector<TwoStateBetaBin> &fg_distro,
                    std::vector<TwoStateBetaBin> &bg_distro) const;

  void
  PosteriorScores(const std::vector<std::vector<std::pair<double, double> > > &values,
                  const std::vector<size_t> &reset_points,
                  const double p_fb, const double p_bf,
                  const std::vector<TwoStateBetaBin> &fg_distro,
                  const std::vector<TwoStateBetaBin> &bg_distro,
                  const bool fg_class,
                  std::vector<double> &llr_scores) const;

  double
  PosteriorDecoding(const std::vector<std::vector<std::pair<double, double> > > &values,
                    const std::vector<size_t> &reset_points,
                    const double p_fb, const double p_bf,
                    const std::vector<TwoStateBetaBin> &fg_distro,
                    const std::vector<TwoStateBetaBin> &bg_distro,
                    std::vector<bool> &classes,
                    std::vector<double> &llr_scores) const;


  double tolerance;
  size_t max_iterations;
  bool VERBOSE;
};

#endif
