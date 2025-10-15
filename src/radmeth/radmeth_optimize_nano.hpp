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

#ifndef RADMETH_OPTIMIZE_NANO_HPP
#define RADMETH_OPTIMIZE_NANO_HPP

#include <vector>

struct RegressionNano;

void
fit_regression_model_nano(RegressionNano &r, std::vector<double> &p_estimates,
                          double &dispersion_estimate);

#endif  // RADMETH_OPTIMIZE_NANO_HPP
