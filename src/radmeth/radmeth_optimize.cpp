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
#include "radmeth_summary_stats.hpp"

#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>

#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <vector>

#include <iostream>

[[nodiscard]] static inline double
logistic(const double x) {
  return 1.0 / (1.0 / std::exp(x) + 1.0);
}

static inline void
logistic(const std::vector<double> &x, std::vector<double> &y) {
  y.resize(std::size(x));
  auto y_itr = std::begin(y);
  for (auto x_itr = std::cbegin(x); x_itr != std::cend(x); ++x_itr)
    *y_itr++ = logistic(*x_itr);
}

static void
set_all_disp_params(const Regression &reg, gsl_vector *disp_params,
                    const double shared_disp) {
  for (auto i = 0u; i < reg.n_disp_params(); ++i)
    gsl_vector_set(disp_params, reg.get_disp_param_idx(i), shared_disp);
}

template <typename T>
[[nodiscard]] static double
get_p(const std::vector<T> &covariates, const gsl_vector_const_view &p_params) {
  // ADS: careful about size of data here
  const auto c = covariates.data();
  return logistic(std::inner_product(c, c + std::size(covariates),
                                     p_params.vector.data, 0.0));
}

template <typename T>
[[nodiscard]] static double
get_param_estimate(const std::vector<T> &covariates,
                   const gsl_vector_const_view &coefficients) {
  const auto c = covariates.data();
  return logistic(std::inner_product(c, c + std::size(covariates),
                                     coefficients.vector.data, 0.0));
}

[[nodiscard]] static double
log_likelihood(Regression &reg, const gsl_vector_const_view &p_coeffs,
               const gsl_vector_const_view &disp_coeffs) {
  // const auto n_groups = reg.n_groups();
  // const auto &groups = reg.design.groups;

  std::vector<double> phi_v(reg.n_disp_params());
  for (std::size_t g_idx = 0; g_idx < reg.n_disp_groups(); ++g_idx) {
    phi_v[g_idx] = get_param_estimate(reg.get_disp_group(g_idx), disp_coeffs);
    std::cerr << "log_likelihood\tphi_v[" << g_idx << "]=" << phi_v[g_idx]
              << std::endl;
  }

  std::vector<double> p_v(reg.n_p_params());
  for (std::size_t g_idx = 0; g_idx < reg.n_p_groups(); ++g_idx) {
    p_v[g_idx] = get_param_estimate(reg.get_p_group(g_idx), p_coeffs);
    std::cerr << "log_likelihood\tp_v[" << g_idx << "]=" << p_v[g_idx]
              << std::endl;
  }

  // cache_log1p_factors(reg, phis);  // ADS: precompute log1p(phi * (k - 1.0))

  const auto &cumul = reg.cumul;
  // std::vector<double> p_params(reg.n_p_params());

  double log_lik = 0.0;
  for (std::size_t g_idx = 0; g_idx < reg.n_groups(); ++g_idx) {
    const auto disp_idx = reg.get_disp_idx(g_idx);
    const auto phi = phi_v[disp_idx];
    const auto one_minus_phi = 1.0 - phi;

    const auto p_idx = reg.get_p_idx(g_idx);
    const auto p = p_v[p_idx];
    const auto one_minus_p = 1.0 - p;

    const auto term1 = one_minus_phi * p;
    const auto &cumul_y = cumul[g_idx].m_counts;
    for (std::size_t k = 0; k < std::size(cumul_y); ++k)
      log_lik += cumul_y[k] * std::log(term1 + phi * k);

    const auto term2 = one_minus_phi * one_minus_p;
    const auto &cumul_d = cumul[g_idx].d_counts;
    for (std::size_t k = 0; k < std::size(cumul_d); ++k)
      log_lik += cumul_d[k] * std::log(term2 + phi * k);

    // const auto &log1p_fact_v = reg.disp_cache[phi_idx];
    const auto &cumul_n = cumul[g_idx].r_counts;
    for (std::size_t k = 0; k < std::size(cumul_n); ++k)
      log_lik -= cumul_n[k] * log1p(phi * (k - 1.0));  // log1p_fact_v[k];
  }
  std::cerr << log_lik << std::endl;
  exit(0);
  return log_lik;
}

static void
gradient(Regression &reg, const gsl_vector_const_view &p_coeffs,
         const gsl_vector_const_view &disp_coeffs, gsl_vector_view &output) {

  std::vector<double> phi_v(reg.n_disp_params());
  for (std::size_t g_idx = 0; g_idx < reg.n_disp_groups(); ++g_idx)
    phi_v[g_idx] = get_param_estimate(reg.get_disp_group(g_idx), disp_coeffs);

  std::vector<double> p_v(reg.n_p_params());
  for (std::size_t g_idx = 0; g_idx < reg.n_p_groups(); ++g_idx)
    p_v[g_idx] = get_param_estimate(reg.get_p_group(g_idx), p_coeffs);

  // init output to zero for all factors
  gsl_vector_set_all(&output.vector, 0.0);
  auto &derivs = output.vector.data;

  // std::vector<double> disp_params(reg.n_disp_params());
  // get_disp_params(reg, params, disp_params);
  // std::vector<double> phis;
  // logistic(disp_params, phis);

  // const auto n_groups = reg.n_groups();
  // const auto &groups = reg.design.groups;
  const auto &cumul = reg.cumul;

  // auto &p_v = reg.p_v;  // ADS: reusing scratch space
  // for (auto g_idx = 0u; g_idx < n_groups; ++g_idx)
  //   p_v[g_idx] = get_p(groups[g_idx], params);

  // cache_dispersion_effects(reg, phis);

  for (std::size_t g_idx = 0; g_idx < reg.n_p_groups(); ++g_idx) {
    const auto disp_idx = reg.get_disp_idx(g_idx);
    const auto phi = phi_v[disp_idx];
    const auto one_minus_phi = 1.0 - phi;

    const auto p_idx = reg.get_p_idx(g_idx);
    const auto p = p_v[p_idx];
    const auto one_minus_p = 1.0 - p;

    double deriv = 0.0;

    const auto denom_term1 = one_minus_phi * p;
    const auto &cumul_y = cumul[g_idx].m_counts;
    const auto y_lim = std::size(cumul_y);
    for (auto k = 0u; k < y_lim; ++k)
      deriv += cumul_y[k] / (denom_term1 + phi * k);

    const auto denom_term2 = one_minus_phi * one_minus_p;
    const auto &cumul_d = cumul[g_idx].d_counts;
    const auto d_lim = std::size(cumul_d);
    for (auto k = 0u; k < d_lim; ++k)
      deriv -= cumul_d[k] / (denom_term2 + phi * k);

    // include the current g_idx group's contribution to each derivative of a
    // 'p' parameter.
    const auto &g = reg.design_p.groups[p_idx];
    const auto denom_term1_one_minus_p = denom_term1 * one_minus_p;
    for (auto p_idx = 0u; p_idx < reg.n_p_params(); ++p_idx) {
      const auto level = g[p_idx];
      derivs[p_idx] += deriv * (denom_term1_one_minus_p * level);
    }
  }

  for (std::size_t g_idx = 0; g_idx < reg.n_disp_groups(); ++g_idx) {
    const auto disp_idx = reg.get_disp_idx(g_idx);
    const auto phi = phi_v[disp_idx];
    const auto one_minus_phi = 1.0 - phi;

    const auto p_idx = reg.get_p_idx(g_idx);
    const auto p = p_v[p_idx];
    const auto one_minus_p = 1.0 - p;

    double disp_deriv = 0.0;

    const auto denom_term1 = one_minus_phi * p;
    const auto &cumul_y = cumul[g_idx].m_counts;
    const auto y_lim = std::size(cumul_y);
    for (auto k = 0u; k < y_lim; ++k)
      disp_deriv += (k - p) * (cumul_y[k] / (denom_term1 + phi * k));

    const auto denom_term2 = one_minus_phi * one_minus_p;
    const auto &cumul_d = cumul[g_idx].d_counts;
    const auto d_lim = std::size(cumul_d);
    for (auto k = 0u; k < d_lim; ++k)
      disp_deriv += (k - one_minus_p) * (cumul_d[k] / (denom_term2 + phi * k));

    const auto &cumul_n = cumul[g_idx].r_counts;
    const auto &disp_effect = reg.disp_cache[g_idx];  // (k-1)/(1 + phi(k-1))
    const auto n_lim = std::size(cumul_n);
    for (auto k = 0u; k < n_lim; ++k)
      disp_deriv -=
        cumul_n[k] * ((k - 1.0) / (1.0 + phi * (k - 1.0)));  // disp_effect[k];

    gsl_vector_set(&output.vector, disp_idx,
                   disp_deriv * (phi * one_minus_phi));
  }
}

[[nodiscard]] static double
neg_loglik(const gsl_vector *params, void *object) {
  auto reg = static_cast<Regression *>(object);
  const auto n_p = reg->n_p_params();
  const auto n_disp = reg->n_disp_params();

  auto p_view = gsl_vector_const_subvector(params, 0, n_p);
  auto disp_view = gsl_vector_const_subvector(params, n_p, n_disp);

  std::vector<double> phi_v(reg->n_disp_params());
  for (std::size_t g_idx = 0; g_idx < reg->n_disp_groups(); ++g_idx) {
    phi_v[g_idx] = get_param_estimate(reg->get_disp_group(g_idx), disp_view);
    std::cerr << "log_likelihood\tphi_v[" << g_idx << "]=" << phi_v[g_idx]
              << std::endl;
  }

  // const gsl_vector *p_params = &p_view.vector;
  // const gsl_vector *disp_params = &disp_view.vector;

  return -log_likelihood(*reg, p_view, disp_view);
}

static void
neg_gradient(const gsl_vector *params, void *object, gsl_vector *output) {
  auto reg = static_cast<Regression *>(object);
  const auto n_p = reg->n_p_params();
  const auto n_disp = reg->n_disp_params();

  auto p_view = gsl_vector_const_subvector(params, 0, n_p);
  auto disp_view = gsl_vector_const_subvector(params, n_p, n_disp);
  auto output_view = gsl_vector_subvector(output, 0, n_p + n_disp);

  // const gsl_vector *p_params = &p_view.vector;
  // const gsl_vector *disp_params = &disp_view.vector;

  gradient(*reg, p_view, disp_view, output_view);
  gsl_vector_scale(output, -1.0);
}

static void
neg_loglik_and_grad(const gsl_vector *params, void *object, double *loglik_val,
                    gsl_vector *output) {
  auto reg = static_cast<Regression *>(object);
  const auto n_p = reg->n_p_params();
  const auto n_disp = reg->n_disp_params();

  auto p_view = gsl_vector_const_subvector(params, 0, n_p);
  auto disp_view = gsl_vector_const_subvector(params, n_p, n_disp);
  auto output_view = gsl_vector_subvector(output, 0, n_p + n_disp);

  // const gsl_vector *p_params = &p_view.vector;
  // const gsl_vector *disp_params = &disp_view.vector;

  *loglik_val = -log_likelihood(*reg, p_view, disp_view);
  gradient(*reg, p_view, disp_view, output_view);
  gsl_vector_scale(&output_view.vector, -1.0);
}

void
fit_regression_model(Regression &reg, std::vector<double> &p_v,
                     std::vector<double> &phi_v) {
  static constexpr auto init_dispersion_param = -2.5;
  const auto stepsize = Regression::stepsize;
  const auto max_iter = Regression::max_iter;

  const auto n_groups = reg.n_groups();
  get_cumulative(reg.get_group_id(), reg.n_groups(), reg.props.mc, reg.cumul);
  set_max_r_count(reg);

  reg.p_v.resize(n_groups);

  const std::size_t n_params = reg.n_params();
  const auto tol =
    std::sqrt(n_params) * reg.n_samples() * Regression::tolerance;
  // clang-format off
  auto loglik_bundle = gsl_multimin_function_fdf{
    &neg_loglik,                // objective function
    &neg_gradient,              // gradient
    &neg_loglik_and_grad,       // combined obj and grad
    n_params,                   // number of model params
    static_cast<void *>(&reg)   // parameters for objective and gradient functions
  };
  // clang-format on

  // set the parameters: zero for "p" parameters and the final one for
  // dispersion using the constant
  auto params = gsl_vector_alloc(n_params);
  gsl_vector_set_all(params, 0.0);
  set_all_disp_params(reg, params, init_dispersion_param);

  // Alternatives:
  // - gsl_multimin_fdfminimizer_conjugate_pr
  // - gsl_multimin_fdfminimizer_conjugate_fr
  // - gsl_multimin_fdfminimizer_vector_bfgs2
  // - gsl_multimin_fdfminimizer_steepest_descent
  const auto minimizer = gsl_multimin_fdfminimizer_vector_bfgs2;
  auto s = gsl_multimin_fdfminimizer_alloc(minimizer, n_params);

  std::cerr << "BEFORE" << std::endl;
  gsl_multimin_fdfminimizer_set(s, &loglik_bundle, params, stepsize, tol);
  std::cerr << "AFTER" << std::endl;
  exit(0);

  int status = 0;
  std::size_t iter = 0;
  do {
    status = gsl_multimin_fdfminimizer_iterate(s);  // one iter and get status
    std::cerr << status << std::endl;
    if (status)
      break;

    const auto param_estimates = gsl_multimin_fdfminimizer_x(s);
    gsl_vector_fprintf(stdout, param_estimates, "%g");

    // check status from gradient
    status = gsl_multimin_test_gradient(s->gradient, tol);
  } while (status == GSL_CONTINUE && ++iter < max_iter);

  /// ADS: the condition below might not always pass even if we are doing
  /// ok. It's not clear how to check for failure.

  // if (status != GSL_SUCCESS)
  //   throw std::runtime_error("failed to fit model parameters");

  const auto param_estimates = gsl_multimin_fdfminimizer_x(s);
  const auto n_p = reg.n_p_params();
  const auto n_disp = reg.n_disp_params();
  auto p_coeffs = gsl_vector_const_subvector(param_estimates, 0, n_p);
  auto disp_coeffs = gsl_vector_const_subvector(param_estimates, n_p, n_disp);

  reg.max_loglik = log_likelihood(reg, p_coeffs, disp_coeffs);
  std::cerr << reg.max_loglik << std::endl;

  phi_v = std::vector<double>(reg.n_disp_params());
  for (std::size_t g_idx = 0; g_idx < reg.n_disp_groups(); ++g_idx) {
    const auto x = get_param_estimate(reg.get_disp_group(g_idx), disp_coeffs);
    phi_v[g_idx] = 1.0 / std::exp(x);
    std::cerr << "phi_v[" << g_idx << "]=" << phi_v[g_idx] << std::endl;
  }

  p_v = std::vector<double>(reg.n_p_params());
  for (std::size_t g_idx = 0; g_idx < reg.n_p_groups(); ++g_idx) {
    p_v[g_idx] = get_param_estimate(reg.get_p_group(g_idx), p_coeffs);
    std::cerr << "p_v[" << g_idx << "]=" << p_v[g_idx] << std::endl;
  }

  gsl_multimin_fdfminimizer_free(s);
  gsl_vector_free(params);
}
