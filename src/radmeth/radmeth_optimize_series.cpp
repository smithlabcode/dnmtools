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

#include "radmeth_optimize_series.hpp"

#include "radmeth_model.hpp"
#include "radmeth_optimize_params.hpp"

#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>

#include <algorithm>
#include <cstring>
#include <numeric>
#include <stdexcept>
#include <vector>

#include <chrono>    /////////////////////
#include <iostream>  /////////////////////

[[nodiscard]] static inline auto
logistic(const double x) -> double {
  return 1.0 / (1.0 / std::exp(x) + 1.0);
}

template <typename T>
[[nodiscard]] static inline auto
get_p(const std::vector<T> &v, const gsl_vector *params) -> double {
  const auto a = v.data();
  return logistic(std::inner_product(a, a + std::size(v), params->data, 0.0));
}

template <typename T>
[[nodiscard]] static inline auto
get_p(const std::vector<T> &v, const std::vector<double> &params) -> double {
  const auto a = v.data();
  return logistic(std::inner_product(a, a + std::size(v), params.data(), 0.0));
}

static inline auto
set_max_r_count(Regression<std::uint32_t> &reg) {
  const auto &cumul = reg.cumul;
  const auto max_itr = std::max_element(
    std::cbegin(cumul), std::cend(cumul), [](const auto &a, const auto &b) {
      return std::size(a.r_counts) < std::size(b.r_counts);
    });
  reg.max_r_count = std::size(max_itr->r_counts);
  // ADS: avoid the realloc that can happen even for resize(smaller_size)
  if (reg.max_r_count > std::size(reg.cache_log1p_factors)) {
    reg.cache_log1p_factors.resize(reg.max_r_count);
    reg.cache_dispersion_effect.resize(reg.max_r_count);
  }
}

static inline auto
get_cached_log1p_factors(Regression<std::uint32_t> &reg, const double phi) {
  const std::size_t max_k = reg.max_r_count;
  auto &cache_log1p_factors = reg.cache_log1p_factors;
  for (std::size_t k = 0; k < max_k; ++k)
    cache_log1p_factors[k] = std::log1p(phi * (k - 1.0));
}

static inline auto
get_cached_dispersion_effect(Regression<std::uint32_t> &reg, const double phi) {
  const std::size_t max_k = reg.max_r_count;
  auto &cache_dispersion_effect = reg.cache_dispersion_effect;
  double j = -1.0;
  const auto lim = std::cbegin(cache_dispersion_effect) + max_k;
  for (auto it = std::begin(cache_dispersion_effect); it != lim; ++it, ++j)
    *it = j / (1.0 + phi * j);
}

[[nodiscard]] static double
log_likelihood(const gsl_vector *params, Regression<std::uint32_t> &reg) {
  const auto phi = logistic(gsl_vector_get(params, reg.design.n_factors()));
  const auto one_minus_phi = 1.0 - phi;

  const auto n_groups = reg.n_groups();
  const auto &groups = reg.design.groups;
  const auto &cumul = reg.cumul;

  // ADS: precompute log1p(phi * (k - 1.0)), reused for each group
  get_cached_log1p_factors(reg, phi);
  const auto &log1p_fact_v = reg.cache_log1p_factors;

  double log_lik = 0.0;
  for (std::size_t g_idx = 0; g_idx < n_groups; ++g_idx) {
    const auto p = get_p(groups[g_idx], params);
    const auto one_minus_p = 1.0 - p;

    const auto term1 = one_minus_phi * p;
    const auto &cumul_y = cumul[g_idx].m_counts;
    std::size_t sz = std::size(cumul_y);
    for (std::size_t k = 0; k < sz; ++k)
      log_lik += cumul_y[k] * std::log(term1 + phi * k);

    const auto term2 = one_minus_phi * one_minus_p;
    const auto &cumul_d = cumul[g_idx].d_counts;
    sz = std::size(cumul_d);
    for (std::size_t k = 0; k < sz; ++k)
      log_lik += cumul_d[k] * std::log(term2 + phi * k);

    const auto &cumul_n = cumul[g_idx].r_counts;
    sz = std::size(cumul_n);
    for (std::size_t k = 0; k < sz; ++k)
      log_lik -= cumul_n[k] * log1p_fact_v[k];
  }
  return log_lik;
}

[[nodiscard]] static double
log_likelihood(const std::vector<double> &params,
               Regression<std::uint32_t> &reg) {
  const auto phi = logistic(params[reg.design.n_factors()]);
  const auto one_minus_phi = 1.0 - phi;

  const auto n_groups = reg.n_groups();
  const auto &groups = reg.design.groups;
  const auto &cumul = reg.cumul;

  // ADS: precompute log1p(phi * (k - 1.0)), reused for each group
  get_cached_log1p_factors(reg, phi);
  const auto &log1p_fact_v = reg.cache_log1p_factors;

  double log_lik = 0.0;
  for (std::size_t g_idx = 0; g_idx < n_groups; ++g_idx) {
    const auto p = get_p(groups[g_idx], params);
    const auto one_minus_p = 1.0 - p;

    const auto term1 = one_minus_phi * p;
    const auto &cumul_y = cumul[g_idx].m_counts;
    std::size_t sz = std::size(cumul_y);
    for (std::size_t k = 0; k < sz; ++k)
      log_lik += cumul_y[k] * std::log(term1 + phi * k);

    const auto term2 = one_minus_phi * one_minus_p;
    const auto &cumul_d = cumul[g_idx].d_counts;
    sz = std::size(cumul_d);
    for (std::size_t k = 0; k < sz; ++k)
      log_lik += cumul_d[k] * std::log(term2 + phi * k);

    const auto &cumul_n = cumul[g_idx].r_counts;
    sz = std::size(cumul_n);
    for (std::size_t k = 0; k < sz; ++k)
      log_lik -= cumul_n[k] * log1p_fact_v[k];
  }
  return -log_lik;
}

static void
gradient(const gsl_vector *params, Regression<std::uint32_t> &reg,
         gsl_vector *output) {
  const auto n_factors = reg.design.n_factors();
  const auto phi = logistic(gsl_vector_get(params, n_factors));
  const auto one_minus_phi = 1.0 - phi;

  const auto n_groups = reg.n_groups();
  const auto &groups = reg.design.groups;
  const auto &cumul = reg.cumul;

  auto &p_v = reg.p_v;
  for (auto g_idx = 0u; g_idx < n_groups; ++g_idx)
    p_v[g_idx] = get_p(groups[g_idx], params);

  get_cached_dispersion_effect(reg, phi);  // (k-1)/(1 + phi(k-1))
  const auto &dispersion_effect = reg.cache_dispersion_effect;

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

static void
gradient(const std::vector<double> &params, Regression<std::uint32_t> &reg,
         std::vector<double> &output) {
  const auto n_factors = reg.design.n_factors();
  const auto phi = logistic(params[n_factors]);
  const auto one_minus_phi = 1.0 - phi;

  const auto n_groups = reg.n_groups();
  const auto &groups = reg.design.groups;
  const auto &cumul = reg.cumul;

  auto &p_v = reg.p_v;
  for (auto g_idx = 0u; g_idx < n_groups; ++g_idx)
    p_v[g_idx] = get_p(groups[g_idx], params);

  get_cached_dispersion_effect(reg, phi);  // (k-1)/(1 + phi(k-1))
  const auto &dispersion_effect = reg.cache_dispersion_effect;

  // init output to zero for all factors
  output.clear();
  output.resize(std::size(params), 0.0);
  auto &data = output;

  double disp_deriv = 0.0;
  for (std::size_t g_idx = 0; g_idx < n_groups; ++g_idx) {
    const auto p = p_v[g_idx];
    const auto one_minus_p = 1.0 - p;

    // if (p == 0.0 || one_minus_p == 0.0 || phi == 0.0 || one_minus_phi == 0.0)
    // {
    //   std::cout << p << '\t' << one_minus_p << std::endl;
    //   std::cout << phi << '\t' << one_minus_phi << std::endl;
    // }

    double deriv = 0.0;

    const auto denom_term1 = one_minus_phi * p;
    const auto &cumul_y = cumul[g_idx].m_counts;
    const auto y_lim = std::size(cumul_y);
    for (auto k = 0u; k < y_lim; ++k) {
      const auto common_factor = cumul_y[k] / (denom_term1 + phi * k);
      deriv += common_factor;
      disp_deriv += (k - p) * common_factor;
      // if (!std::isfinite(disp_deriv)) {
      //   // clang-format off
      //   std::cout << "common_factor: " << common_factor << '\n'
      //             << "cumul_y[k]: " << cumul_y[k] << '\n'
      //             << "denom_term1: " << denom_term1 << '\n'
      //             << "phi: " << phi << '\n'
      //             << "k: " << k << '\n'
      //             << "p: " << p << '\n'
      //             << "one_minus_phi: " << one_minus_phi << '\n'
      //             << "phi: " << phi << '\n'
      //             << "params[n_factors]: " << params[n_factors] << '\n'
      //             << "(denom_term1 + phi * k): " << (denom_term1 + phi * k)
      //             << '\n'
      //             << std::endl;
      //   throw std::runtime_error("bad disp_deriv 1");
      //   // clang-format on
      // }
    }

    const auto denom_term2 = one_minus_phi * one_minus_p;
    const auto &cumul_d = cumul[g_idx].d_counts;
    const auto d_lim = std::size(cumul_d);
    for (auto k = 0u; k < d_lim; ++k) {
      const auto common_factor = cumul_d[k] / (denom_term2 + phi * k);
      deriv -= common_factor;
      disp_deriv += (k - one_minus_p) * common_factor;
      // if (!std::isfinite(disp_deriv))
      //   throw std::runtime_error("bad disp_deriv 2");
    }

    const auto &cumul_n = cumul[g_idx].r_counts;
    const auto n_lim = std::size(cumul_n);
    for (auto k = 0u; k < n_lim; ++k) {
      disp_deriv -= cumul_n[k] * dispersion_effect[k];
      // if (!std::isfinite(disp_deriv))
      //   throw std::runtime_error("bad disp_deriv 3");
    }

    const auto &g = groups[g_idx];
    const auto denom_term1_one_minus_p = denom_term1 * one_minus_p;
    for (auto fact_idx = 0u; fact_idx < n_factors; ++fact_idx) {
      const auto level = g[fact_idx];
      if (level == 0)
        continue;
      data[fact_idx] += deriv * (denom_term1_one_minus_p * level);
    }
  }
  // if (!std::isfinite(disp_deriv))
  //   throw std::runtime_error("bad disp_deriv");
  output[n_factors] = disp_deriv * (phi * one_minus_phi);
  // for (auto i = 0u; i <= n_factors; ++i) {
  //   if (!std::isfinite(output[i]))
  //     throw std::runtime_error("bad2");
  // }
  for (auto i = 0u; i <= n_factors; ++i) {
    output[i] *= -1.0;
    // if (!std::isfinite(output[i]))
    //   throw std::runtime_error("bad");
  }
  // std::cout << "grad=" << disp_deriv * (phi * one_minus_phi) << std::endl;
}

[[nodiscard]] static double
neg_loglik(const gsl_vector *params, void *object) {
  auto reg = static_cast<Regression<std::uint32_t> *>(object);
  return -log_likelihood(params, *reg);
}

static void
neg_gradient(const gsl_vector *params, void *object, gsl_vector *output) {
  auto reg = static_cast<Regression<std::uint32_t> *>(object);
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
               const std::uint32_t n_groups,
               const std::vector<mcounts<std::uint32_t>> &mc,
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
  // call 3 times: m_counts, r_counts, d_counts
  comp_cumul(
    [](const mcounts<std::uint32_t> &m) { return m.n_meth; },
    [](cumul_counts &c) -> std::vector<std::uint32_t> & { return c.m_counts; });

  comp_cumul(
    [](const mcounts<std::uint32_t> &m) { return m.n_reads; },
    [](cumul_counts &c) -> std::vector<std::uint32_t> & { return c.r_counts; });

  comp_cumul(
    [](const mcounts<std::uint32_t> &m) { return m.n_reads - m.n_meth; },
    [](cumul_counts &c) -> std::vector<std::uint32_t> & { return c.d_counts; });
}

static inline auto
conjugate_pr(std::vector<double> &params, Regression<std::uint32_t> &reg,
             const double alpha) -> int {
  const auto n = std::size(params);
  const auto max_iter = radmeth_optimize_params::max_iter;
  const auto tol = radmeth_optimize_params::tolerance;

  std::vector<double> g(n);
  std::vector<double> g_prev(n);
  gradient(params, reg, g);

  auto d = g;
  std::transform(std::begin(d), std::end(d), std::begin(d),
                 [](const auto x) { return -1.0 * x; });

  // const auto ll = log_likelihood(params, reg);
  // std::cout << "ll=" << ll << std::endl;

  double gdot{};
  for (const double gi : g)
    gdot += gi * gi;
  double denom = gdot;
  double gnorm = std::sqrt(gdot);

  std::size_t iter{};
  for (iter = 0; iter < max_iter; ++iter) {
    // std::memcpy(g_prev.data(), g.data(), sizeof(double) * n);
    // update parameters based on directions and step size
    for (std::size_t i = 0; i < n; ++i)
      params[i] += alpha * d[i];

    gradient(params, reg, g);

    // get numerator and for norm of gradient
    gdot = 0.0;
    for (const double gi : g)
      gdot += gi * gi;
    gnorm = std::sqrt(gdot);

    // stop if norm of gradient is small enough
    if (gnorm < tol)
      break;

    // // get denom
    // double denom = 0.0;
    // for (const double gpi : g_prev)
    //   denom += gpi * gpi;

    // update conjugate directions
    const double beta = std::max(0.0, gdot / denom);
    for (auto i = 0u; i < n; ++i)
      d[i] = -g[i] + beta * d[i];

    denom = gdot;
  }
  return iter;
  // if (gnorm < tol) {
  //   std::cout << "Converged in " << iter << " iterations" << "\t" << gnorm
  //             << '\t' << tol << std::endl;
  // }
  // else {
  //   std::cout << "Did not converge in " << max_iter << " iterations\n";
  // }
}

void
fit_regression_model(Regression<std::uint32_t> &r,
                     std::vector<double> &p_estimates,
                     double &dispersion_estimate) {
  static constexpr auto init_dispersion_param = -2.5;
  const auto stepsize = radmeth_optimize_params::stepsize;
  const auto max_iter = radmeth_optimize_params::max_iter;

  const auto n_groups = r.n_groups();
  get_cumulative(r.design.group_id, n_groups, r.mc, r.cumul);
  set_max_r_count(r);

  r.p_v.resize(n_groups);

  const std::size_t n_params = r.n_params();
  const auto tol =
    std::sqrt(n_params) * r.n_samples() * radmeth_optimize_params::tolerance;
  // clang-format off
  auto loglik_bundle = gsl_multimin_function_fdf{
    &neg_loglik,              // objective function
    &neg_gradient,            // gradient
    &neg_loglik_and_grad,     // combined obj and grad
    n_params,                 // number of model params
    static_cast<void *>(&r)   // params for objective and gradients
  };
  // clang-format on

  // std::vector<double> x(n_params, 0.0);
  // // for (auto i = 0u; i < n_params; ++i)
  // //   x[i] = gsl_vector_get(param_estimates, i);
  // x.back() = -2.5;
  // bool all_finite = false;
  // auto start_time_a = std::chrono::steady_clock::now();
  // auto ii = 0;
  // for (ii = 0; ii < 4 && !all_finite; ++ii) {
  //   double alpha = std::pow(10.0, -1.0 * (ii + 2));
  //   std::cout << "alpha=" << alpha << std::endl;
  //   std::vector<double> x_out(x);
  //   conjugate_pr(x_out, r, alpha);
  //   all_finite = true;
  //   for (auto i = 0u; i < std::size(x_out) && all_finite; ++i)
  //     all_finite = std::isfinite(x_out[i]);
  //   if (all_finite)
  //     x = x_out;
  // }
  // // if (!all_finite) {
  // //   alpha = 1e-03;
  // //   x_out = x;
  // //   ii = conjugate_pr(x, r, alpha);
  // //   all_finite = true;
  // //   for (auto i = 0u; i < std::size(x) && all_finite; ++i)
  // //     all_finite = std::isfinite(x[i]);
  // //   if (!all_finite) {
  // //     alpha = 1e-04;
  // //     x_out = x;
  // //     ii = conjugate_pr(x, r, alpha);
  // //     all_finite = true;
  // //     for (auto i = 0u; i < std::size(x) && all_finite; ++i)
  // //       all_finite = std::isfinite(x[i]);
  // //     if (!all_finite) {
  // //       alpha = 1e-05;
  // //       x_out = x;
  // //       ii = conjugate_pr(x, r, alpha);
  // //       all_finite = true;
  // //       for (auto i = 0u; i < std::size(x) && all_finite; ++i)
  // //         all_finite = std::isfinite(x[i]);
  // //       if (!all_finite) {
  // //         alpha = 1e-06;
  // //         x_out = x;
  // //         ii = conjugate_pr(x, r, alpha);
  // //         all_finite = true;
  // //         for (auto i = 0u; i < std::size(x) && all_finite; ++i)
  // //           all_finite = std::isfinite(x[i]);
  // //       }
  // //     }
  // //   }
  // // }
  // auto stop_time_a = std::chrono::steady_clock::now();
  // // std::cout << "[total time: " << (stop_time_a - start_time_a).count() <<
  // // "]\t"
  // //           << ii << std::endl;

  // parameters: 0 for "p" params and final one is a constant for dispersion
  auto params = gsl_vector_alloc(n_params);
  gsl_vector_set_all(params, 0.0);
  gsl_vector_set(params, n_params - 1, init_dispersion_param);

  // Alternatives:
  // - gsl_multimin_fdfminimizer_conjugate_pr
  // - gsl_multimin_fdfminimizer_conjugate_fr
  // - gsl_multimin_fdfminimizer_vector_bfgs2
  // - gsl_multimin_fdfminimizer_steepest_descent

  const auto minimizer = gsl_multimin_fdfminimizer_conjugate_fr;
  auto s = gsl_multimin_fdfminimizer_alloc(minimizer, n_params);
  gsl_multimin_fdfminimizer_set(s, &loglik_bundle, params, stepsize, tol);

  // gsl_multimin_function minex_func;

  // minex_func.n = n_params;
  // minex_func.f = &neg_loglik;
  // minex_func.params = static_cast<void *>(&r);

  // /* Set initial step sizes to 1 */
  // auto ss = gsl_vector_alloc(n_params);
  // gsl_vector_set_all(ss, 1.0);

  // const auto minimizer = gsl_multimin_fminimizer_nmsimplex2;
  // auto s = gsl_multimin_fminimizer_alloc(minimizer, n_params);
  // gsl_multimin_fminimizer_set(s, &minex_func, params, ss);

  auto start_time_b = std::chrono::steady_clock::now();
  int status = 0;
  std::size_t iter = 0;

  double size{};
  do {
    status = gsl_multimin_fdfminimizer_iterate(s);  // one iter and get
    // status
    // status = gsl_multimin_fminimizer_iterate(s);  // one iter and get status
    if (status)
      break;
    // check status from gradient

    status = gsl_multimin_test_gradient(s->gradient, tol);

    // size = gsl_multimin_fminimizer_size(s);
    // status = gsl_multimin_test_size(size, 1e-2);

  } while (status == GSL_CONTINUE && ++iter < max_iter);
  // ADS: can't use (status != GSL_SUCCESS)
  auto stop_time_b = std::chrono::steady_clock::now();
  // std::cout << "[total time: "
  //           << (stop_time_b - start_time_b).count()
  //           // -
  //           //              (stop_time_a - start_time_a).count()
  //           << "]\t" << iter << std::endl;

  const auto param_estimates = gsl_multimin_fdfminimizer_x(s);
  // const auto param_estimates = gsl_multimin_fminimizer_x(s);
  // for (auto i = 0u; i < n_params; ++i)
  //   std::cout << "x[i]=" << gsl_vector_get(param_estimates, i) << std::endl;

  // gsl_vector_fprintf(stdout, param_estimates, "%g");

  const auto &groups = r.design.groups;
  p_estimates.clear();
  for (auto g_idx = 0u; g_idx < n_groups; ++g_idx)
    p_estimates.push_back(get_p(groups[g_idx], param_estimates));
  const auto disp_param = gsl_vector_get(param_estimates, n_params - 1);
  dispersion_estimate = 1.0 / std::exp(disp_param);

  r.max_loglik = log_likelihood(param_estimates, r);

  gsl_multimin_fdfminimizer_free(s);
  // gsl_multimin_fminimizer_free(s);
  gsl_vector_free(params);
}
