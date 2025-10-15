/* Copyright (C) 2013-2025 Andrew D Smith
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

#include "radmeth_optimize_nano.hpp"
#include "radmeth_model_nano.hpp"

#include <gsl/gsl_multimin.h>
#include <gsl/gsl_sf_psi.h>  // for gsl_sf_psi (digamma)
#include <gsl/gsl_vector.h>

#include <algorithm>
#include <array>
#include <numeric>
#include <stdexcept>
#include <vector>

/* Coefficients for the Chebyschev polynomial for the digamma function in the
   range 0-1 */
// clang-format off
static constexpr std::array<double, 23> psi1_cs {
  -0.038057080835217922, // == -0.019028540417608961*2
   0.491415393029387130,
  -0.056815747821244730,
   0.008357821225914313,
  -0.001333232857994342,
   0.000220313287069308,
  -0.000037040238178456,
   0.000006283793654854,
  -0.000001071263908506,
   0.000000183128394654,
  -0.000000031353509361,
   0.000000005372808776,
  -0.000000000921168141,
   0.000000000157981265,
  -0.000000000027098646,
   0.000000000004648722,
  -0.000000000000797527,
   0.000000000000136827,
  -0.000000000000023475,
   0.000000000000004027,
  -0.000000000000000691,
   0.000000000000000118,
  -0.000000000000000020
};
// clang-format on

/* Alternate set of coefficients for the Chebyschev polynomial for the digamma
   function */
// clang-format off
static constexpr std::array<double, 16> psi2_cs = {
  -0.0204749044678185, // == -0.01023745223390925*2
  -0.0101801271534859,
   0.0000559718725387,
  -0.0000012917176570,
   0.0000000572858606,
  -0.0000000038213539,
   0.0000000003397434,
  -0.0000000000374838,
   0.0000000000048990,
  -0.0000000000007344,
   0.0000000000001233,
  -0.0000000000000228,
   0.0000000000000045,
  -0.0000000000000009,
   0.0000000000000002,
  -0.0000000000000000 ,
};
// clang-format on

template <std::size_t T>
[[nodiscard]] static inline double
chebyschev(const std::array<double, T> &coeffs, std::uint32_t order,
           const double y) {
  const auto y2 = 2.0 * y;
  double d = 0.0;
  double dd = 0.0;
  for (auto j = order; j >= 1; j--) {
    const auto temp = d;
    d = y2 * d - dd + coeffs[j];
    dd = temp;
  }
  return y * d - dd + 0.5 * coeffs[0];
}

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

static void
get_cache_lgamma(const std::vector<std::uint8_t> &group,
                 const gsl_vector *params, const double phi,
                 vars_cache &cache) {
  const auto p = get_p(group, params);
  const auto a = p * phi;
  const auto b = (1.0 - p) * phi;
  cache.p = p;
  cache.a = a;
  cache.b = b;
  cache.lgamma_a = std::lgamma(a);
  cache.lgamma_b = std::lgamma(b);
  cache.lgamma_a_b = std::lgamma(a + b);
}

[[nodiscard]] static double
log_likelihood_nano(const gsl_vector *params, RegressionNano &reg) {
  const auto phi = 1.0 / std::exp(gsl_vector_get(params, reg.n_factors()));

  const auto n_groups = reg.n_groups();
  const auto &groups = reg.design.groups;
  auto &cache = reg.cache;
  for (auto g_idx = 0u; g_idx < n_groups; ++g_idx)
    get_cache_lgamma(groups[g_idx], params, phi, cache[g_idx]);

  const auto &group_id = reg.design.group_id;
  const auto &mc = reg.mc;
  double ll = 0.0;
  const auto n_samples = reg.n_samples();
  for (auto i = 0u; i < n_samples; ++i) {
    const auto y = mc[i].n_meth;
    const auto n = mc[i].n_reads;

    const auto &c = cache[group_id[i]];
    const auto a = c.a;
    const auto b = c.b;

    // clang-format off
    ll += ((std::lgamma(y + a) - c.lgamma_a) +
           (std::lgamma(n - y + b) - c.lgamma_b) +
           (c.lgamma_a_b - std::lgamma(n + a + b)));
    // clang-format on
  }

  return ll;
}

// digamma for x non-negative
static double
digamma(const double y) {
  static constexpr auto psi_order = 7;   // max=22;
  static constexpr auto apsi_order = 7;  // max=15;
  if (y >= 2.0) {
    const auto t = 8.0 / (y * y) - 1.0;
    return std::log(y) - 0.5 / y + chebyschev(psi2_cs, apsi_order, t);
  }
  if (y < 1.0)
    return -1.0 / y + chebyschev(psi1_cs, psi_order, 2.0 * y - 1.0);
  return chebyschev(psi1_cs, psi_order, 2.0 * (y - 1.0) - 1.0);
}

static void
get_cache_digamma(const std::vector<std::uint8_t> &group,
                  const gsl_vector *params, const double phi,
                  vars_cache &cache) {
  const auto p = get_p(group, params);
  const auto a = p * phi;
  const auto b = (1.0 - p) * phi;
  cache.p = p;
  cache.a = a;
  cache.b = b;
  cache.digamma_a = digamma(a);
  cache.digamma_b = digamma(b);
  cache.digamma_a_b = digamma(a + b);
}

static void
gradient_nano(const gsl_vector *params, RegressionNano &reg,
              gsl_vector *output) {
  const auto n_factors = reg.n_factors();
  const auto phi = 1.0 / std::exp(gsl_vector_get(params, n_factors));

  const auto &mc = reg.mc;
  const auto &matrix = reg.design.matrix;

  const auto n_groups = reg.n_groups();
  const auto &groups = reg.design.groups;
  auto &cache = reg.cache;  // ADS: reusing scratch space
  for (auto g_idx = 0u; g_idx < n_groups; ++g_idx)
    get_cache_digamma(groups[g_idx], params, phi, cache[g_idx]);

  gsl_vector_set_all(output, 0.0);  // init output to zero for all factors
  auto &data = output->data;

  const auto &group_id = reg.design.group_id;
  double grad_phi = 0.0;
  const auto n_samples = reg.n_samples();
  for (auto i = 0u; i < n_samples; ++i) {
    const auto y = mc[i].n_meth;
    const auto n = mc[i].n_reads;

    const auto &c = cache[group_id[i]];
    const auto p = c.p;
    const auto a = c.a;
    const auto b = c.b;

    const auto digamma_a = c.digamma_a;
    const auto digamma_b = c.digamma_b;
    const auto digamma_y_a = digamma(y + a);
    const auto digamma_n_y_b = digamma(n - y + b);

    // grad wrt p
    const auto dlogl_dp =
      phi * (digamma_y_a - digamma_a - digamma_n_y_b + digamma_b);
    const auto dp_delta = p * (1.0 - p);           // chain rule: d p / d param
    const auto dlogl_delta = dlogl_dp * dp_delta;  // grad wrt params

    auto matrix_itr = std::cbegin(matrix[i]);
    const auto data_end = data + n_factors;
    for (auto data_itr = data; data_itr != data_end; ++data_itr)
      *data_itr += (*matrix_itr++) * dlogl_delta;

    const auto digamma_delta = c.digamma_a_b - digamma(n + a + b);
    const auto dphi_term1 = p * (digamma_y_a - digamma_a + digamma_delta);
    const auto dphi_term2 =
      (1.0 - p) * (digamma_n_y_b - digamma_b + digamma_delta);

    grad_phi += dphi_term1 + dphi_term2;
  }

  const auto grad_theta = -grad_phi * phi;
  gsl_vector_set(output, n_factors, grad_theta);
}

[[nodiscard]] static double
neg_loglik_nano(const gsl_vector *params, void *object) {
  auto reg = static_cast<RegressionNano *>(object);
  return -log_likelihood_nano(params, *reg);
}

static void
neg_gradient_nano(const gsl_vector *params, void *object, gsl_vector *output) {
  auto reg = static_cast<RegressionNano *>(object);
  gradient_nano(params, *reg, output);
  gsl_vector_scale(output, -1.0);
}

static void
neg_loglik_and_grad_nano(const gsl_vector *params, void *object,
                         double *loglik_val, gsl_vector *d_loglik_val) {
  *loglik_val = neg_loglik_nano(params, object);
  neg_gradient_nano(params, object, d_loglik_val);
}

void
fit_regression_model_nano(RegressionNano &r, std::vector<double> &p_estimates,
                          double &dispersion_estimate) {
  static constexpr auto init_dispersion_param = -2.5;
  const auto stepsize = RegressionNano::stepsize;
  const auto max_iter = RegressionNano::max_iter;
  const auto n_groups = r.n_groups();
  r.cache.resize(n_groups);  // make sure scratch space is allocated

  const std::size_t n_params = r.n_params();
  const auto tol =
    std::sqrt(n_params) * r.n_samples() * RegressionNano::tolerance;
  // clang-format off
  auto loglik_bundle = gsl_multimin_function_fdf{
    &neg_loglik_nano,              // objective function
    &neg_gradient_nano,            // gradient
    &neg_loglik_and_grad_nano,     // combined obj and grad
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

  r.max_loglik = log_likelihood_nano(param_estimates, r);

  gsl_multimin_fdfminimizer_free(s);
  gsl_vector_free(params);
}
