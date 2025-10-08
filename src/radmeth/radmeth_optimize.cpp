/* Copyright (C) 2013-2025 Andrew D Smith
 *
 * Author: Andrew D. Smith
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
#include <gsl/gsl_sf_psi.h>  // for gsl_sf_psi (digamma)
#include <gsl/gsl_vector.h>

#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <vector>

[[nodiscard]] static inline double
logistic(const double x) {
  return 1.0 / (1.0 / std::exp(x) + 1.0);
}

template <typename T>
[[nodiscard]] static double
get_p(const std::vector<T> &v, const gsl_vector *params) {
  const auto a = v.data();
  return logistic(std::inner_product(a, a + std::size(v), params->data, 0.0));
}

[[nodiscard]] static double
log_likelihood(const gsl_vector *params, Regression &reg) {
  const auto phi = 1.0 / std::exp(gsl_vector_get(params, reg.n_factors()));

  const auto &mc = reg.props.mc;
  const auto &matrix = reg.design.matrix;

  double ll = 0.0;
  const auto n_samples = reg.n_samples();
  for (auto i = 0u; i < n_samples; ++i) {
    const auto p = get_p(matrix[i], params);
    const auto y = mc[i].n_meth;
    const auto n = mc[i].n_reads;

    const auto a = p * phi;
    const auto b = (1.0 - p) * phi;

    // clang-format off
    ll += ((std::lgamma(y + a) - std::lgamma(a)) +
           (std::lgamma(n - y + b) - std::lgamma(b)) +
           (std::lgamma(a + b) - std::lgamma(n + a + b)));
    // clang-format on
  }

  return ll;
}

[[nodiscard]] static inline double
digamma(const double x) {
  return gsl_sf_psi(x);
}

// gradient contribution for one obs wrt p
[[nodiscard]] static double
grad_loglik_p(const double y, const std::uint32_t n, const double p,
              const double phi) {
  const auto a = p * phi;
  const auto b = (1.0 - p) * phi;
  return phi * (digamma(y + a) - digamma(a) - digamma(n - y + b) + digamma(b));
}

[[nodiscard]] static double
grad_loglik_phi_single_obs(const double y, const std::uint32_t n,
                           const double p, const double phi) {
  const auto a = p * phi;
  const auto b = (1.0 - p) * phi;

  const auto digamma_aplusb = digamma(a + b);
  const auto digamma_n_aplusb = digamma(n + a + b);
  const auto dg_delta = digamma_aplusb - digamma_n_aplusb;

  const auto term1 = p * (digamma(y + a) - digamma(a) + dg_delta);
  const auto term2 = (1.0 - p) * (digamma(n - y + b) - digamma(b) + dg_delta);

  return term1 + term2;
}

static void
gradient(const gsl_vector *params, Regression &reg, gsl_vector *output) {
  const auto n_factors = reg.n_factors();
  const auto phi = 1.0 / std::exp(gsl_vector_get(params, n_factors));

  const auto &mc = reg.props.mc;
  const auto &matrix = reg.design.matrix;

  gsl_vector_set_all(output, 0.0);  // init output to zero for all factors
  auto &data = output->data;

  double grad_phi = 0.0;
  const auto n_samples = reg.n_samples();
  for (auto i = 0u; i < n_samples; ++i) {
    const auto &matrix_i = matrix[i];
    const auto p = get_p(matrix_i, params);
    const auto y = mc[i].n_meth;
    const auto n = mc[i].n_reads;

    const auto dlogl_dp = grad_loglik_p(y, n, p, phi);  // grad wrt p
    const auto dp_deta = p * (1.0 - p);          // chain rule: d p / d param
    const auto dlogl_deta = dlogl_dp * dp_deta;  // grad wrt params

    for (auto fact_idx = 0u; fact_idx < n_factors; ++fact_idx) {
      const auto level = matrix_i[fact_idx];
      if (level == 0)
        continue;
      data[fact_idx] += dlogl_deta;
    }
    grad_phi += grad_loglik_phi_single_obs(y, n, p, phi);
  }

  const auto grad_theta = -grad_phi * phi;
  gsl_vector_set(output, n_factors, grad_theta);
}

[[nodiscard]] static double
neg_loglik(const gsl_vector *params, void *object) {
  auto reg = static_cast<Regression *>(object);
  return -log_likelihood(params, *reg);
}

static void
neg_gradient(const gsl_vector *params, void *object, gsl_vector *output) {
  auto reg = static_cast<Regression *>(object);
  gradient(params, *reg, output);
  gsl_vector_scale(output, -1.0);
}

static void
neg_loglik_and_grad(const gsl_vector *params, void *object, double *loglik_val,
                    gsl_vector *d_loglik_val) {
  *loglik_val = neg_loglik(params, object);
  neg_gradient(params, object, d_loglik_val);
}

void
fit_regression_model(Regression &r, std::vector<double> &p_estimates,
                     double &dispersion_estimate) {
  static constexpr auto init_dispersion_param = -2.5;
  const auto stepsize = Regression::stepsize;
  const auto max_iter = Regression::max_iter;
  const auto n_groups = r.n_groups();

  r.p_v.resize(n_groups);

  const std::size_t n_params = r.n_params();
  const auto tol = std::sqrt(n_params) * r.n_samples() * Regression::tolerance;
  // clang-format off
  auto loglik_bundle = gsl_multimin_function_fdf{
    &neg_loglik,              // objective function
    &neg_gradient,            // gradient
    &neg_loglik_and_grad,     // combined obj and grad
    n_params,                 // number of model params
    static_cast<void *>(&r)   // parameters for objective and gradient functions
  };
  // clang-format on

  // set the parameters: zero for "p" parameters and the final one for
  // dispersion using the constant
  auto params = gsl_vector_alloc(n_params);
  gsl_vector_set_all(params, 0.0);
  gsl_vector_set(params, n_params - 1, init_dispersion_param);

  // Alternatives:
  // - gsl_multimin_fdfminimizer_conjugate_pr
  // - gsl_multimin_fdfminimizer_conjugate_fr
  // - gsl_multimin_fdfminimizer_vector_bfgs2
  // - gsl_multimin_fdfminimizer_steepest_descent
  const auto minimizer = gsl_multimin_fdfminimizer_conjugate_pr;
  auto s = gsl_multimin_fdfminimizer_alloc(minimizer, n_params);

  gsl_multimin_fdfminimizer_set(s, &loglik_bundle, params, stepsize, tol);

  int status = 0;
  std::size_t iter = 0;
  do {
    status = gsl_multimin_fdfminimizer_iterate(s);  // one iter and get status
    if (status)
      break;
    // check status from gradient
    status = gsl_multimin_test_gradient(s->gradient, tol);
  } while (status == GSL_CONTINUE && ++iter < max_iter);

  /// ADS: the condition below might not always pass even if we are doing
  /// ok. It's not clear how to check for failure.

  // if (status != GSL_SUCCESS)
  //   throw std::runtime_error("failed to fit model parameters");

  const auto param_estimates = gsl_multimin_fdfminimizer_x(s);

  const auto &groups = r.design.groups;
  p_estimates.clear();
  for (auto g_idx = 0u; g_idx < n_groups; ++g_idx)
    p_estimates.push_back(get_p(groups[g_idx], param_estimates));
  const auto disp_param = gsl_vector_get(param_estimates, n_params - 1);
  dispersion_estimate = 1.0 / std::exp(disp_param);

  r.max_loglik = log_likelihood(param_estimates, r);

  gsl_multimin_fdfminimizer_free(s);
  gsl_vector_free(params);
}
