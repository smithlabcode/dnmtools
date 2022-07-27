/*
  Copyright (C) 2022 Andrew D Smith

  Authors: Andrew D. Smith and Benjamin E. Decato

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

#ifndef EM_DTN
#define EM_DTN

#include <cmath>
#include <iomanip>
#include <numeric>
#include <iostream>
#include <sstream>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_gamma.h>
#include <utility>
#include <string>
#include <vector>

/** Emission distributions for methylation should be modeled either as
 * Beta or Beta Binomial. Since they will be used simultaneously, it is
 * helpful to have an abstraction so that we can put them in the same
 * container.
 */
class EmissionDistribution
{
  public:
    EmissionDistribution();
    virtual ~EmissionDistribution();
    EmissionDistribution(const double a, const double b);
    EmissionDistribution(const std::string &str);
    virtual double operator()(const std::pair<double, double> &val) const = 0;
    virtual double log_likelihood(const std::pair<double, double> &val) const = 0;
    std::string tostring() const;
    double getalpha() { return alpha; };
    double getbeta() { return beta; };
    void fit(const std::vector<double> &vals_a,
             const std::vector<double> &vals_b,
             const std::vector<double> &p);

  protected:
    double sign(const double x);
    double invpsi(const double tolerance, const double x);
    double movement(const double curr, const double prev);
    double alpha;
    double beta;
    double lnbeta_helper;

    const double tolerance = 1e-10;
};

class Beta : public EmissionDistribution
{
  public:
    Beta();
    Beta(const double a, const double b);
    Beta(const std::string &str);
    double operator()(const std::pair<double, double> &val) const;
    double log_likelihood(const std::pair<double, double> &val) const;
};

class BetaBinomial : public EmissionDistribution
{
  public:
    BetaBinomial();
    BetaBinomial(const double a, const double b);
    BetaBinomial(const std::string &str);
    double operator()(const std::pair<double, double> &val) const;
    double log_likelihood(const std::pair<double, double> &val) const;
};

#endif
