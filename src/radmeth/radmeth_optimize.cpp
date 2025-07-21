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

[[nodiscard]] static double
pi(const std::vector<double> &v, const gsl_vector *params) {
  // ADS: this function doesn't have a very helpful name
  const auto a = v.data();
  const double dot = std::inner_product(a, a + std::size(v), params->data, 0.0);
  // const double p = -1.0 / std::expm1(-dot_prod); ?
  return std::exp(dot) / (1.0 + std::exp(dot));
}

[[nodiscard]] static double
neg_loglik(const gsl_vector *params, void *object) {
  Regression *reg = (Regression *)(object);

  // ADS: the dispersion parameter phi is the last element of
  // parameter vector
  const double disp_param = gsl_vector_get(params, reg->design.n_factors());
  const double phi = std::exp(disp_param) / (1.0 + std::exp(disp_param));
  const double one_minus_phi = 1.0 - phi;

  const auto &mc = reg->props.mc;

  double log_lik = 0;
  const std::size_t n_samples = reg->design.n_samples();
  for (std::size_t col_idx = 0; col_idx < n_samples; ++col_idx) {
    const double n = mc[col_idx].n_reads;
    const double y = mc[col_idx].n_meth;
    const double p = pi(reg->design.matrix[col_idx], params);
    const double one_minus_p = 1.0 - p;

    const double term1 = one_minus_phi * p;
    for (int k = 0; k < y; ++k)
      log_lik += std::log(term1 + phi * k);

    const double term2 = one_minus_phi * one_minus_p;
    for (int k = 0; k < n - y; ++k)
      log_lik += std::log(term2 + phi * k);

    for (int k = 0; k < n; ++k)
      log_lik -= std::log(1.0 + phi * (k - 1.0));
  }
  return -log_lik;
}

static void
neg_gradient(const gsl_vector *params, void *object, gsl_vector *output) {
  Regression *reg = (Regression *)(object);
  const std::size_t n_samples = reg->design.n_samples();
  const std::size_t n_factors = reg->design.n_factors();

  std::vector<double> p_v(n_samples, 0.0);
  for (auto col_idx = 0u; col_idx < n_samples; ++col_idx)
    p_v[col_idx] = pi(reg->design.matrix[col_idx], params);

  const double disp_param = gsl_vector_get(params, n_factors);

  const auto &tmatrix = reg->design.tmatrix;
  const double phi = std::exp(disp_param) / (1.0 + std::exp(disp_param));
  // const double phi = -1.0 / std::expm1(-disp_param);
  const double one_minus_phi = 1.0 - phi;

  const auto &mc = reg->props.mc;

  for (std::size_t f = 0; f < n_factors; ++f) {
    double deriv = 0;
    const auto &vec = tmatrix[f];
    for (std::size_t col_idx = 0; col_idx < n_samples; ++col_idx) {
      const double n = mc[col_idx].n_reads;
      const double y = mc[col_idx].n_meth;
      const double p = p_v[col_idx];
      const double one_minus_p = 1.0 - p;

      const double denom_term1 = one_minus_phi * p;
      const double factor = denom_term1 * one_minus_p * vec[col_idx];
      if (factor == 0)
        continue;

      double term = 0;

      const auto accumulate_term = [&](int k, const int lim,
                                       const double denom_base) {
        for (; k < lim; ++k)
          term += 1.0 / (denom_base + phi * k);
      };

      const auto accumulate_subtract_term = [&](int k, const int lim,
                                                const double denom_base) {
        for (; k < lim; ++k)
          term -= 1.0 / (denom_base + phi * k);
      };

      accumulate_term(0, y, denom_term1);
      accumulate_subtract_term(0, n - y, one_minus_phi * one_minus_p);

      deriv += term * factor;
    }
    gsl_vector_set(output, f, deriv);
  }

  double deriv = 0;
  for (std::size_t col_idx = 0; col_idx < n_samples; ++col_idx) {
    const double n = mc[col_idx].n_reads;
    const double y = mc[col_idx].n_meth;
    const double p = p_v[col_idx];
    const double one_minus_p = 1.0 - p;

    double term = 0;
    const auto accumulate_term = [&](int k, const int lim,
                                     const double num_shift,
                                     const double denom_base) {
      for (; k < lim; ++k)
        term += (k - num_shift) / (denom_base + phi * k);
    };
    const auto accumulate_subtract_term = [&](int k, const int lim) {
      for (; k < lim; ++k)
        term -= (k - 1.0) / (1.0 + phi * (k - 1.0));
    };

    accumulate_term(0, y, p, one_minus_phi * p);
    accumulate_term(0, n - y, one_minus_p, one_minus_phi * one_minus_p);
    accumulate_subtract_term(0, n);

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
  const auto tolerance = Regression::tolerance;
  const auto stepsize = Regression::stepsize;
  const auto max_iter = Regression::max_iter;

  // one more than the number of factors
  const std::size_t n_params = r.n_factors() + 1;
  if (params_init.empty()) {
    params_init.resize(n_params, 0.0);
    params_init.back() = -2.5;
  }

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
  // It it reasonable to reduce the number of iterations to 500?
  // ADS: 700 vs. 500? what's the difference?

  r.max_loglik = -1.0 * neg_loglik(s->x, &r);

  gsl_multimin_fdfminimizer_free(s);
  gsl_vector_free(params);

  return status == GSL_SUCCESS;
}
