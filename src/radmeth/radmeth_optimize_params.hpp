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

#ifndef RADMETH_OPTIMIZE_PARAMS_HPP
#define RADMETH_OPTIMIZE_PARAMS_HPP

#include <cstdint>

namespace radmeth_optimize_params {
inline double tolerance = 1e-4;
inline double stepsize = 0.01;
inline std::uint32_t max_iter = 250;
};  // namespace radmeth_optimize_params

#endif  // RADMETH_OPTIMIZE_PARAMS_HPP
