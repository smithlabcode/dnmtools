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

#include "radmeth_optimize.hpp"

#include "radmeth_model.hpp"

#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>

#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <vector>

template <typename T>
[[nodiscard]] static double
pi(const std::vector<T> &v, const gsl_vector *params) {
  // ADS: this function doesn't have a very helpful name
  const auto a = v.data();
  const double dot = std::inner_product(a, a + std::size(v), params->data, 0.0);
  // const double p = std::exp(dot) / (1.0 + std::exp(dot));
  return 1.0 / (1.0 / std::exp(dot) + 1.0);
}

[[nodiscard]] static double
neg_loglik(const gsl_vector *params, void *object) {
  Regression *reg = (Regression *)(object);

  // ADS: the dispersion parameter phi is the last element of
  // parameter vector
  const double disp_param = gsl_vector_get(params, reg->design.n_factors());
  // const double phi = std::exp(disp_param) / (1.0 + std::exp(disp_param));
  const double phi = 1.0 / (1.0 / std::exp(disp_param) + 1.0);
  const double one_minus_phi = 1.0 - phi;

  const auto &mc = reg->props.mc;
  const auto &mat = reg->design.matrix;

  double log_lik = 0;
  const std::size_t n_samples = reg->design.n_samples();
  for (std::size_t col_idx = 0; col_idx < n_samples; ++col_idx) {
    const auto n = mc[col_idx].n_reads;
    const auto y = mc[col_idx].n_meth;
    const double p = pi(mat[col_idx], params);
    const double one_minus_p = 1.0 - p;

    const double term1 = one_minus_phi * p;
    for (auto k = 0u; k < y; ++k)
      log_lik += std::log(term1 + phi * k);

    const double term2 = one_minus_phi * one_minus_p;
    for (auto k = 0u; k < n - y; ++k)
      log_lik += std::log(term2 + phi * k);

    for (auto k = 0u; k < n; ++k)
      // log_lik -= std::log(1.0 + phi * (k - 1.0));
      log_lik -= std::log1p(phi * (k - 1.0));
  }
  return -log_lik;
}

static void
neg_gradient(const gsl_vector *params, void *object, gsl_vector *output) {
  Regression *reg = (Regression *)(object);
  const std::size_t n_samples = reg->design.n_samples();
  const std::size_t n_factors = reg->design.n_factors();

  /// ADS: using scratch space held in reg instead of allocating here
  // std::vector<double> p_v(n_samples, 0.0);
  auto &p_v = reg->p_v;
  for (auto col_idx = 0u; col_idx < n_samples; ++col_idx)
    p_v[col_idx] = pi(reg->design.matrix[col_idx], params);

  const double disp_param = gsl_vector_get(params, n_factors);

  const auto &tmat = reg->design.tmatrix;
  // const double phi = std::exp(disp_param) / (1.0 + std::exp(disp_param));
  const double phi = 1.0 / (1.0 / std::exp(disp_param) + 1.0);
  const double one_minus_phi = 1.0 - phi;

  const auto &mc = reg->props.mc;

  for (std::size_t factor_idx = 0; factor_idx < n_factors; ++factor_idx) {
    double deriv = 0;
    const auto &vec = tmat[factor_idx];
    for (std::size_t col_idx = 0; col_idx < n_samples; ++col_idx) {
      const auto n = mc[col_idx].n_reads;
      const auto y = mc[col_idx].n_meth;
      const double p = p_v[col_idx];
      const double one_minus_p = 1.0 - p;

      const double denom_term1 = one_minus_phi * p;
      const double factor = denom_term1 * one_minus_p * vec[col_idx];
      if (factor == 0)
        continue;

      double term = 0;

      const auto accum_term = [&](auto k, const auto lim,
                                  const double denom_base) {
        for (; k < lim; ++k)
          term += 1.0 / (denom_base + phi * k);
      };

      const auto accum_subtract_term = [&](auto k, const auto lim,
                                           const double denom_base) {
        for (; k < lim; ++k)
          term -= 1.0 / (denom_base + phi * k);
      };

      accum_term(0u, y, denom_term1);
      accum_subtract_term(0u, n - y, one_minus_phi * one_minus_p);

      deriv += term * factor;
    }
    gsl_vector_set(output, factor_idx, deriv);
  }

  double deriv = 0;
  for (std::size_t col_idx = 0; col_idx < n_samples; ++col_idx) {
    const auto n = mc[col_idx].n_reads;
    const auto y = mc[col_idx].n_meth;
    const double p = p_v[col_idx];
    const double one_minus_p = 1.0 - p;

    double term = 0;
    const auto accum_term = [&](auto k, const auto lim, const double num_shift,
                                const double denom_base) {
      for (; k < lim; ++k)
        term += (k - num_shift) / (denom_base + phi * k);
    };
    const auto accum_subtract_term = [&](auto k, const auto lim) {
      for (; k < lim; ++k)
        term -= (k - 1.0) / (1.0 + phi * (k - 1.0));
    };

    accum_term(0u, y, p, one_minus_phi * p);
    accum_term(0u, n - y, one_minus_p, one_minus_phi * one_minus_p);
    accum_subtract_term(0u, n);

    deriv += term * phi * one_minus_phi;
  }
  gsl_vector_set(output, n_factors, deriv);
  gsl_vector_scale(output, -1.0);
}

static void
neg_loglik_and_grad(const gsl_vector *params, void *object, double *loglik_val,
                    gsl_vector *d_loglik_val) {
  *loglik_val = neg_loglik(params, object);
  neg_gradient(params, object, d_loglik_val);
}

bool
fit_regression_model(Regression &r, std::vector<double> &params_init) {
  const auto stepsize = Regression::stepsize;
  const auto max_iter = Regression::max_iter;

  // one more than the number of factors
  const std::size_t n_params = r.n_factors() + 1;
  if (params_init.empty()) {
    params_init.resize(n_params, 0.0);
    params_init.back() = -2.5;
  }

  r.p_v.resize(r.n_samples());

  const double tolerance =
    std::sqrt(n_params) * r.n_samples() * Regression::tolerance;

  if (params_init.size() != n_params)
    throw std::runtime_error("Wrong number of initial parameters.");

  // clang-format off
  gsl_multimin_function_fdf loglik_bundle = {
    &neg_loglik,              // objective function
    &neg_gradient,            // gradient
    &neg_loglik_and_grad,     // combined obj and grad
    n_params,                 // number of model params
    static_cast<void *>(&r)   // parameters for objective and gradient functions
  };
  // clang-format on

  // Alternatives:
  // - gsl_multimin_fdfminimizer_conjugate_pr
  // - gsl_multimin_fdfminimizer_conjugate_fr
  // - gsl_multimin_fdfminimizer_vector_bfgs2
  // - gsl_multimin_fdfminimizer_steepest_descent
  const gsl_multimin_fdfminimizer_type *minimizer =
    gsl_multimin_fdfminimizer_conjugate_fr;

  gsl_multimin_fdfminimizer *s =
    gsl_multimin_fdfminimizer_alloc(minimizer, n_params);

  gsl_vector *params = gsl_vector_alloc(n_params);
  for (std::size_t i = 0; i < n_params; ++i)
    gsl_vector_set(params, i, params_init[i]);

  gsl_multimin_fdfminimizer_set(s, &loglik_bundle, params, stepsize, tolerance);

  int status = 0;
  std::size_t iter = 0;

  do {
    status = gsl_multimin_fdfminimizer_iterate(s);  // one iter and get status
    if (status)
      break;
    // check status from gradient
    status = gsl_multimin_test_gradient(s->gradient, tolerance);
  } while (status == GSL_CONTINUE && ++iter < max_iter);
  // ADS: What is a reasonable number of iterations?

  r.max_loglik = -1.0 * neg_loglik(s->x, &r);

  gsl_multimin_fdfminimizer_free(s);
  gsl_vector_free(params);

  return status == GSL_SUCCESS;
}
