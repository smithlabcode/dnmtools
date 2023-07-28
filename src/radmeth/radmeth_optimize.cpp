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

#include "radmeth_model.hpp"

#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>

#include <stdexcept>
#include <vector>

using std::vector;
using std::runtime_error;

static double
pi(Regression *reg, size_t sample, const gsl_vector *params) {
  // ADS: this function doesn't have a very helpful name
  double dot_prod = 0;
  for (size_t fact = 0; fact < reg->design.n_factors(); ++fact)
    dot_prod += reg->design.matrix[sample][fact]*gsl_vector_get(params, fact);

  double p = exp(dot_prod)/(1 + exp(dot_prod));

  return p;
}


static double
neg_loglik(const gsl_vector *params, void *object) {
  Regression *reg = (Regression *)(object);
  const size_t n_params = reg->design.factor_names.size() + 1;

  double log_lik = 0;

  // ADS: the dispersion parameter phi is the last element of
  // parameter vector
  const double disp_param = gsl_vector_get(params, n_params - 1);
  const double phi = exp(disp_param)/(1.0 + exp(disp_param));

  for(size_t s = 0; s < reg->design.n_samples(); ++s) {
    const double n_s = reg->props.total[s];
    const double y_s = reg->props.meth[s];
    const double p_s = pi(reg, s, params);

    for (int k = 0; k < y_s; ++k)
      log_lik += log((1 - phi)*p_s + phi*k);

    for (int k = 0; k < n_s - y_s; ++k)
      log_lik += log((1 - phi)*(1 - p_s) + phi*k);

    for (int k = 0; k < n_s; ++k)
      log_lik -= log(1.0 + phi*(k - 1));
  }

  return -log_lik;
}


static void
neg_gradient(const gsl_vector *params, void *object,
             gsl_vector *output) {

  Regression *reg = (Regression *)(object);
  const size_t n_params = reg->design.n_factors() + 1;

  const double disp_param = gsl_vector_get(params, n_params - 1);

  const double phi = exp(disp_param)/(1 + exp(disp_param));

  for(size_t f = 0; f < n_params; ++f) {

    double deriv = 0;

    for(size_t s = 0; s < reg->design.n_samples(); ++s) {
      int n_s = reg->props.total[s];
      int y_s = reg->props.meth[s];
      double p_s = pi(reg, s, params);

      double term = 0;

      //a parameter linked to p
      if(f < reg->design.factor_names.size()) {
        double factor = (1 - phi)*p_s*(1 - p_s)*reg->design.matrix[s][f];
        if (factor == 0) continue;

        for(int k = 0; k < y_s; ++k)
          term += 1/((1 - phi)*p_s + phi*k);

        for(int k = 0; k < n_s - y_s; ++k)
          term -= 1/((1 - phi)*(1 - p_s) + phi*k);

        deriv += term*factor;
      } else { // the parameter linked to phi
        for(int k = 0; k < y_s; ++k)
          term += (k - p_s)/((1 - phi)*p_s + phi*k);

        for(int k = 0; k < n_s - y_s; ++k)
          term += (k - (1 - p_s))/((1 - phi)*(1 - p_s) + phi*k);

        for(int k = 0; k < n_s; ++k) {
          term -= (k - 1)/(1 + phi*(k - 1));
        }

        deriv += term * phi * (1 - phi);
      }
    }

    gsl_vector_set(output, f, deriv);
  }

  gsl_vector_scale(output, -1.0);
}

static void
neg_loglik_and_grad(const gsl_vector *params,
                    void *object,
                    double *loglik_val,
                    gsl_vector *d_loglik_val) {

  *loglik_val = neg_loglik(params, object);
  neg_gradient(params, object, d_loglik_val);
}


bool
fit_regression_model(Regression &r, vector<double> &initial_params) {

  // one more than the number of factors
  const size_t n_params = r.n_factors() + 1;
  if (initial_params.empty()) {
    initial_params.resize(n_params, 0.0);
    initial_params.back() = -2.5;
  }

  if (initial_params.size() != n_params)
    throw runtime_error("Wrong number of initial parameters.");

  int status = 0;

  size_t iter = 0;

  {
    gsl_multimin_function_fdf loglik_bundle;
    loglik_bundle.f = &neg_loglik;
    loglik_bundle.df = &neg_gradient;
    loglik_bundle.fdf = &neg_loglik_and_grad;
    loglik_bundle.n = n_params;
    loglik_bundle.params = (void *)&r;

    gsl_vector *params = gsl_vector_alloc(n_params);

    for (size_t param = 0; param < initial_params.size(); ++param)
      gsl_vector_set(params, param, initial_params[param]);

    //can also try gsl_multimin_fdfminimizer_conjugate_pr;
    const gsl_multimin_fdfminimizer_type *T = gsl_multimin_fdfminimizer_conjugate_fr;

    gsl_multimin_fdfminimizer *s = gsl_multimin_fdfminimizer_alloc (T, n_params);

    // ADS: 0.001 is the step size and 1e-4 is the tolerance
    // get the minimizer
    gsl_multimin_fdfminimizer_set(s, &loglik_bundle, params, 0.001, 1e-4);

    do {
      iter++;
      // do one iteration and get the status
      status = gsl_multimin_fdfminimizer_iterate(s);
      if (status)
        break;
      // check the status based on gradient
      status = gsl_multimin_test_gradient (s->gradient, 1e-4);
    }
    while (status == GSL_CONTINUE && iter < 700);
    // It it reasonable to reduce the number of iterations to 500?
    // ADS: 700 vs. 500? what's the difference?

    r.max_loglik = (-1)*neg_loglik(s->x, &r);

    gsl_multimin_fdfminimizer_free(s);
    gsl_vector_free(params);
  }

  return status == GSL_SUCCESS;
}
