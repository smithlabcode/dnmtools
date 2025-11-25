/*
  Copyright (C) 2011-2022 University of Southern California
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

#include "BetaBin.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <utility>

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>

using std::cerr;
using std::max;
using std::min;
using std::pair;
using std::setprecision;
using std::setw;
using std::string;
using std::vector;

//////////////////////////////////////////////
//////       struct betabin             //////
//////////////////////////////////////////////

const double betabin::tolerance = 1e-10;

betabin::betabin() : alpha(1), beta(1), lnbeta_helper(gsl_sf_lnbeta(1, 1)) {}

betabin::betabin(const double a, const double b) :
  alpha(a), beta(b), lnbeta_helper(gsl_sf_lnbeta(a, b)) {}

betabin::betabin(const string &str) {
  std::istringstream iss(str, std::istringstream::in);
  string name;
  iss >> name >> alpha >> beta;
  if (name != "betabin" || alpha < 0 || beta < 0) {
    cerr << "betabin::betabin: "
         << "bad string representation of betabin distribution: " << str
         << '\n';
    throw "bad string representation of betabin distribution";
  }
  lnbeta_helper = gsl_sf_lnbeta(alpha, beta);
}

string
betabin::tostring() const {
  std::ostringstream os;
  os << "betabin " << setprecision(4) << alpha << " " << setprecision(4)
     << beta;
  return os.str();
}

double
betabin::operator()(const pair<double, double> &val) const {
  const std::uint32_t x = static_cast<std::uint32_t>(val.first);
  const std::uint32_t n = static_cast<std::uint32_t>(val.first + val.second);
  return gsl_sf_lnchoose(n, x) + gsl_sf_lnbeta(alpha + x, beta + val.second) -
         lnbeta_helper;
}

double
betabin::log_likelihood(const pair<double, double> &val) const {
  const std::uint32_t x = static_cast<std::uint32_t>(val.first);
  const std::uint32_t n = static_cast<std::uint32_t>(val.first + val.second);
  return gsl_sf_lnchoose(n, x) + gsl_sf_lnbeta(alpha + x, beta + val.second) -
         lnbeta_helper;
}

double
betabin::sign(const double x) {
  return (x >= 0) ? 1.0 : -1.0;
}

double
betabin::invpsi(const double tolerance, const double x) {
  double L = 1.0;
  double Y = std::exp(x);
  while (L > tolerance) {
    Y += L * sign(x - gsl_sf_psi(Y));
    L /= 2.0;  // NOLINT(*-avoid-magic-numbers)
  }
  return Y;
}

double
betabin::movement(const double curr, const double prev) {
  return std::abs(curr - prev) / std::max(std::fabs(curr), std::fabs(prev));
}

void
betabin::fit(const vector<double> &vals_a, const vector<double> &vals_b,
             const vector<double> &p) {
  static constexpr auto initial_param_vals = 0.01;
  const double p_total = std::accumulate(p.begin(), p.end(), 0.0);
  const double alpha_rhs =
    std::inner_product(std::cbegin(vals_a), std::cend(vals_a), std::cbegin(p),
                       0.0) /
    p_total;
  const double beta_rhs =
    std::inner_product(std::cbegin(vals_b), std::cend(vals_b), std::cbegin(p),
                       0.0) /
    p_total;
  double prev_alpha = 0.0, prev_beta = 0.0;
  alpha = beta = initial_param_vals;

  while (movement(alpha, prev_alpha) > tolerance &&
         movement(beta, prev_beta) > tolerance) {
    prev_alpha = alpha;
    prev_beta = beta;
    alpha = invpsi(tolerance, gsl_sf_psi(prev_alpha + prev_beta) + alpha_rhs);
    beta = invpsi(tolerance, gsl_sf_psi(prev_alpha + prev_beta) + beta_rhs);
  }
  lnbeta_helper = gsl_sf_lnbeta(alpha, beta);
}
