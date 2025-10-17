/* Copyright (C) 2025 Andrew D. Smith
 *
 * Author: Andrew D. Smith
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 */

#include "radmeth_optimize.hpp"
#include "radmeth_model.hpp"

#include <gsl/gsl_multimin.h>
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

static inline auto
set_max_r_count(Regression &reg) {
  const auto &cumul = reg.cumul;
  const auto max_itr = std::max_element(
    std::cbegin(cumul), std::cend(cumul), [](const auto &a, const auto &b) {
      return std::size(a.r_counts) < std::size(b.r_counts);
    });
  reg.max_r_count = std::size(max_itr->r_counts);
  // ADS: avoid the realloc that can happen even for resize(smaller_size)
  if (reg.max_r_count > std::size(reg.cache))
    reg.cache.resize(reg.max_r_count);
}

static inline auto
cache_log1p_factors(Regression &reg, const double phi) {
  const std::size_t max_k = reg.max_r_count;
  auto &cache = reg.cache;
  for (std::size_t k = 0; k < max_k; ++k)
    cache[k] = std::log1p(phi * (k - 1.0));
}

static inline auto
cache_dispersion_effect(Regression &reg, const double phi) {
  const std::size_t max_k = reg.max_r_count;
  auto &cache = reg.cache;
  for (std::size_t k = 0; k < max_k; ++k)
    cache[k] = (k - 1.0) / (1.0 + phi * (k - 1.0));
}

[[nodiscard]] static double
log_likelihood(const gsl_vector *params, Regression &reg) {
  const auto phi = logistic(gsl_vector_get(params, reg.design.n_factors()));
  const auto one_minus_phi = 1.0 - phi;

  const auto n_groups = reg.n_groups();
  const auto &groups = reg.design.groups;
  const auto &cumul = reg.cumul;

  // ADS: precompute the log1p(phi * (k - 1.0)) values, which are reused for
  // each group.
  cache_log1p_factors(reg, phi);
  const auto &log1p_fact_v = reg.cache;

  double log_lik = 0.0;
  for (std::size_t g_idx = 0; g_idx < n_groups; ++g_idx) {
    const auto p = get_p(groups[g_idx], params);
    const auto one_minus_p = 1.0 - p;

    const auto term1 = one_minus_phi * p;
    const auto &cumul_y = cumul[g_idx].m_counts;
    for (std::size_t k = 0; k < std::size(cumul_y); ++k)
      log_lik += cumul_y[k] * std::log(term1 + phi * k);

    const auto term2 = one_minus_phi * one_minus_p;
    const auto &cumul_d = cumul[g_idx].d_counts;
    for (std::size_t k = 0; k < std::size(cumul_d); ++k)
      log_lik += cumul_d[k] * std::log(term2 + phi * k);

    const auto &cumul_n = cumul[g_idx].r_counts;
    for (std::size_t k = 0; k < std::size(cumul_n); ++k)
      log_lik -= cumul_n[k] * log1p_fact_v[k];
  }
  return log_lik;
}

static void
gradient(const gsl_vector *params, Regression &reg, gsl_vector *output) {
  const auto n_factors = reg.design.n_factors();
  const auto phi = logistic(gsl_vector_get(params, n_factors));
  const auto one_minus_phi = 1.0 - phi;

  const auto n_groups = reg.n_groups();
  const auto &groups = reg.design.groups;
  const auto &cumul = reg.cumul;

  auto &p_v = reg.p_v;  // ADS: reusing scratch space
  for (auto g_idx = 0u; g_idx < n_groups; ++g_idx)
    p_v[g_idx] = get_p(groups[g_idx], params);

  cache_dispersion_effect(reg, phi);
  const auto &dispersion_effect = reg.cache;  // (k-1)/(1 + phi(k-1))

  // init output to zero for all factors
  gsl_vector_set_all(output, 0.0);
  const auto &data = output->data;

  double disp_deriv = 0.0;
  for (std::size_t g_idx = 0; g_idx < n_groups; ++g_idx) {
    const auto p = p_v[g_idx];
    const auto one_minus_p = 1.0 - p;

    double deriv = 0.0;

    const auto denom_term1 = one_minus_phi * p;
    const auto &cumul_y = cumul[g_idx].m_counts;
    const auto y_lim = std::size(cumul_y);
    for (auto k = 0u; k < y_lim; ++k) {
      const auto common_factor = cumul_y[k] / (denom_term1 + phi * k);
      deriv += common_factor;
      disp_deriv += (k - p) * common_factor;
    }

    const auto denom_term2 = one_minus_phi * one_minus_p;
    const auto &cumul_d = cumul[g_idx].d_counts;
    const auto d_lim = std::size(cumul_d);
    for (auto k = 0u; k < d_lim; ++k) {
      const auto common_factor = cumul_d[k] / (denom_term2 + phi * k);
      deriv -= common_factor;
      disp_deriv += (k - one_minus_p) * common_factor;
    }

    const auto &cumul_n = cumul[g_idx].r_counts;
    const auto n_lim = std::size(cumul_n);
    for (auto k = 0u; k < n_lim; ++k)
      disp_deriv -= cumul_n[k] * dispersion_effect[k];

    const auto &g = groups[g_idx];
    const auto denom_term1_one_minus_p = denom_term1 * one_minus_p;
    for (auto fact_idx = 0u; fact_idx < n_factors; ++fact_idx) {
      const auto level = g[fact_idx];
      if (level == 0)
        continue;
      data[fact_idx] += deriv * (denom_term1_one_minus_p * level);
    }
  }
  gsl_vector_set(output, n_factors, disp_deriv * (phi * one_minus_phi));
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

static void
get_cumulative(const std::vector<std::uint32_t> &group_id,
               const std::uint32_t n_groups, const std::vector<mcounts> &mc,
               std::vector<cumul_counts> &cumul) {
  const auto n_cols = std::size(mc);
  cumul.clear();
  cumul.resize(n_groups);

  const auto comp_cumul = [&](auto get_value, auto get_vector) {
    // phase 1: determine max value for each group
    for (auto g_idx = 0u; g_idx < n_groups; ++g_idx) {
      std::uint32_t max_v{};
      for (auto c_idx = 0u; c_idx < n_cols; ++c_idx) {
        if (group_id[c_idx] == g_idx) {
          const auto val = get_value(mc[c_idx]);
          if (val > max_v)
            max_v = val;
        }
      }
      get_vector(cumul[g_idx]).resize(max_v, 0);
    }

    // phase 2: fill cumulative counts
    for (auto c_idx = 0u; c_idx < n_cols; ++c_idx) {
      const auto g_idx = group_id[c_idx];
      const auto val = get_value(mc[c_idx]);
      auto &vec = get_vector(cumul[g_idx]);
      for (auto i = 0u; i < val; ++i)
        ++vec[i];
    }
  };
  // call the lambda 3 times for m_counts, r_counts, d_counts
  comp_cumul(
    [](const mcounts &m) { return m.n_meth; },
    [](cumul_counts &c) -> std::vector<std::uint32_t> & { return c.m_counts; });

  comp_cumul(
    [](const mcounts &m) { return m.n_reads; },
    [](cumul_counts &c) -> std::vector<std::uint32_t> & { return c.r_counts; });

  comp_cumul(
    [](const mcounts &m) { return m.n_reads - m.n_meth; },
    [](cumul_counts &c) -> std::vector<std::uint32_t> & { return c.d_counts; });
}

void
fit_regression_model(Regression &r, std::vector<double> &p_estimates,
                     double &dispersion_estimate) {
  static constexpr auto init_dispersion_param = -2.5;
  const auto stepsize = Regression::stepsize;
  const auto max_iter = Regression::max_iter;

  const auto n_groups = r.n_groups();
  get_cumulative(r.design.group_id, n_groups, r.mc, r.cumul);
  set_max_r_count(r);

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
