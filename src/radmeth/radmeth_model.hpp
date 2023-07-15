/* Copyright (C) 2013-2023 University of Southern California and
 *                         Egor Dolzhenko
 *                         Andrew D Smith
 *                         Guilherme Sena
 *
 * Authors: Andrew D. Smith and Egor Dolzhenko and Guilherme Sena
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

#ifndef RADMETH_MODEL_HPP
#define RADMETH_MODEL_HPP

#include <string>
#include <vector>

struct Design {
  std::vector<std::string> factor_names;
  std::vector<std::string> sample_names;
  std::vector<std::vector<double> > matrix;
  size_t n_factors() const {return factor_names.size();}
  size_t n_samples() const {return sample_names.size();}
};

struct SiteProportions {
  std::string chrom;
  size_t position;
  std::string strand;
  std::string context;
  std::vector<size_t> total;
  std::vector<size_t> meth;
};

struct Regression {
  Design design;
  SiteProportions props;
  double max_loglik;
  size_t n_factors() const {return design.n_factors();}
  size_t n_samples() const {return design.n_samples();}
};

#endif
