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

#ifndef BETABIN_HPP
#define BETABIN_HPP

#include <utility>
#include <string>
#include <vector>

struct betabin {
  betabin();
  betabin(const double a, const double b);
  betabin(const std::string &str);
  double operator()(const std::pair<double, double> &val) const;
  double log_likelihood(const std::pair<double, double> &val) const;
  double sign(const double x);
  double invpsi(const double tolerance, const double x);
  double movement(const double curr, const double prev);
  void fit(const std::vector<double> &vals_a,
           const std::vector<double> &vals_b,
           const std::vector<double> &p);
  std::string tostring() const;
  double alpha;
  double beta;
  double lnbeta_helper;

  static const double tolerance;
};

#endif
