/* Copyright (C) 2025 Andrew D Smith
 *
 * Author: Andrew D. Smith
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 */

#ifndef RADMETH_OPTIMIZE_GAMMA_HPP
#define RADMETH_OPTIMIZE_GAMMA_HPP

#include <cstdint>
#include <vector>

template <typename T> struct Regression;

void
fit_regression_model_gamma(Regression<std::uint32_t> &r,
                           std::vector<double> &p_estimates,
                           double &dispersion_estimate);

void
fit_regression_model_gamma(Regression<double> &r,
                           std::vector<double> &p_estimates,
                           double &dispersion_estimate);

void
fit_regression_model_gamma_fdf(Regression<std::uint32_t> &r,
                               std::vector<double> &p_estimates,
                               double &dispersion_estimate);

void
fit_regression_model_gamma_fdf(Regression<double> &r,
                               std::vector<double> &p_estimates,
                               double &dispersion_estimate);

#endif  // RADMETH_OPTIMIZE_GAMMA_HPP
