/* Code adapted from GSL, see copyright below.
 */

/* cdf/inverse_normal.c
 *
 * Copyright (C) 2002 Przemyslaw Sliwa and Jason H. Stover.
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc., 51
 * Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/*
 * Computes the inverse normal cumulative distribution function
 * according to the algorithm shown in
 *
 *      Wichura, M.J. (1988).
 *      Algorithm AS 241: The Percentage Points of the Normal Distribution.
 *      Applied Statistics, 37, 477-484.
 */

// AND ALSO:

/* cdf/gauss.c
 *
 * Copyright (C) 2002, 2004 Jason H. Stover.
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc., 51
 * Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/*
 * Computes the cumulative distribution function for the Gaussian distribution
 * using a rational function approximation.  The computation is for the
 * standard Normal distribution, i.e., mean 0 and standard deviation 1. If you
 * want to compute Pr(X < t) for a Gaussian random variable X with non-zero
 * mean m and standard deviation sd not equal to 1, find gsl_cdf_ugaussian
 * ((t-m)/sd).  This approximation is accurate to at least double
 * precision. The accuracy was verified with a pari-gp script.  The largest
 * error found was about 1.4E-20. The coefficients were derived by Cody.
 *
 * References:
 *
 * W.J. Cody. "Rational Chebyshev Approximations for the Error
 * Function," Mathematics of Computation, v23 n107 1969, 631-637.
 *
 * W. Fraser, J.F Hart. "On the Computation of Rational Approximations
 * to Continuous Functions," Communications of the ACM, v5 1962.
 *
 * W.J. Kennedy Jr., J.E. Gentle. "Statistical Computing." Marcel Dekker. 1980.
 *
 *
 */

#include "dnmtools_gaussinv.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <limits>
#include <numeric>

// ADS: so much magic.

template <std::size_t a_size, std::size_t b_size>
[[nodiscard]] static inline double
rat_eval(const std::array<double, a_size> &a,
         const std::array<double, b_size> &b, const double x) {
  const auto acc = [&](const auto t, const auto y) { return x * t + y; };
  return std::accumulate(std::cbegin(a), std::cend(a), 0.0, acc) /
         std::accumulate(std::cbegin(b), std::cend(b), 0.0, acc);
}

[[nodiscard]] static double
small(const double q) {
  // clang-format off
  static constexpr auto a = std::array{
    2509.0809287301226727,
    33430.575583588128105,
    67265.770927008700853,
    45921.953931549871457,
    13731.693765509461125,
    1971.5909503065514427,
    133.14166789178437745,
    3.387132872796366608,
  };
  static constexpr auto b = std::array{
    5226.495278852854561,
    28729.085735721942674,
    39307.89580009271061,
    21213.794301586595867,
    5394.1960214247511077,
    687.1870074920579083,
    42.313330701600911252,
    1.0,
  };
  // clang-format on

  const double r = 0.180625 - q * q;

  return q * rat_eval(a, b, r);
}

[[nodiscard]] static double
intermediate(const double r) {
  static constexpr auto offset = 1.6;
  // clang-format off
  static constexpr auto a = std::array{
    7.7454501427834140764e-4,
    0.0227238449892691845833,
    0.24178072517745061177,
    1.27045825245236838258,
    3.64784832476320460504,
    5.7694972214606914055,
    4.6303378461565452959,
    1.42343711074968357734,
  };
  static constexpr auto b = std::array{
    1.05075007164441684324e-9,
    5.475938084995344946e-4,
    0.0151986665636164571966,
    0.14810397642748007459,
    0.68976733498510000455,
    1.6763848301838038494,
    2.05319162663775882187,
    1.0,
  };
  // clang-format on

  return rat_eval(a, b, r - offset);
}

[[nodiscard]] static double
tail(const double r) {
  static constexpr auto offset = 0.5;
  // clang-format off
  static constexpr auto a = std::array{
    2.01033439929228813265e-7,
    2.71155556874348757815e-5,
    0.0012426609473880784386,
    0.026532189526576123093,
    0.29656057182850489123,
    1.7848265399172913358,
    5.4637849111641143699,
    6.6579046435011037772,
  };
  static constexpr auto b = std::array{
    2.04426310338993978564e-15,
    1.4215117583164458887e-7,
    1.8463183175100546818e-5,
    7.868691311456132591e-4,
    0.0148753612908506148525,
    0.13692988092273580531,
    0.59983220655588793769,
    1.0,
  };
  // clang-format on

  return rat_eval(a, b, r - offset);
}

[[nodiscard]] double
dnmt_gsl_cdf_ugaussian_Pinv(const double P) {
  static constexpr auto cutoff_for_small = 0.425;
  const auto dP = P - 0.5;
  if (P == 1.0)
    return std::numeric_limits<double>::infinity();
  if (P == 0.0)
    return -std::numeric_limits<double>::infinity();
  if (std::abs(dP) <= cutoff_for_small)
    return small(dP);
  const double pp = P < 0.5 ? P : 1.0 - P;
  const double r = std::sqrt(-std::log(pp));
  const double x = (r <= 5.0) ? intermediate(r) : tail(r);
  return P < 0.5 ? -x : x;  // NOLINT(*-avoid-magic-numbers)
}

[[nodiscard]] double
dnmt_gsl_cdf_ugaussian_Qinv(const double Q) {
  static constexpr auto cutoff_for_small = 0.425;
  const double dQ = Q - 0.5;
  if (Q == 1.0)
    return -std::numeric_limits<double>::infinity();
  if (Q == 0.0)
    return std::numeric_limits<double>::infinity();
  if (std::abs(dQ) <= cutoff_for_small)
    return -small(dQ);
  const double pp = Q < 0.5 ? Q : 1.0 - Q;
  const double r = std::sqrt(-std::log(pp));
  const double x = (r <= 5.0) ? intermediate(r) : tail(r);
  return Q < 0.5 ? x : -x;  // NOLINT(*-avoid-magic-numbers)
}

[[nodiscard]] double
dnmt_gsl_cdf_gaussian_Pinv(const double P, const double sigma) {
  return sigma * dnmt_gsl_cdf_ugaussian_Pinv(P);
}

[[nodiscard]] double
dnmt_gsl_cdf_gaussian_Qinv(const double Q, const double sigma) {
  return sigma * dnmt_gsl_cdf_ugaussian_Qinv(Q);
}

// ADS: from gsl/gsl_math.h:

// clang-format off
static constexpr const double dnmt_M_2_SQRTPI = 1.12837916709551257389615890312;  // 2/sqrt(pi)
static constexpr const double dnmt_M_SQRT1_2 = 0.70710678118654752440084436210;  // sqrt(1/2)
static constexpr const double dnmt_M_1_SQRT2PI = dnmt_M_2_SQRTPI * dnmt_M_SQRT1_2 / 2.0;
static constexpr const double dnmt_M_SQRT2 = 1.41421356237309504880168872421;  // sqrt(2)
static constexpr const double dnmt_SQRT32 = 4.0 * dnmt_M_SQRT2;
// clang-format on

/*
 * IEEE double precision dependent constants.
 *
 * GAUSS_EPSILON: Smallest positive value such that gsl_cdf_gaussian(x) > 0.5.
 * GAUSS_XUPPER: Largest value x such that gsl_cdf_gaussian(x) < 1.0.
 * GAUSS_XLOWER: Smallest value x such that gsl_cdf_gaussian(x) > 0.0.
 */

static constexpr const double dnmt_GSL_DBL_EPSILON = 2.2204460492503131e-16;
static constexpr const double dnmt_GAUSS_EPSILON = dnmt_GSL_DBL_EPSILON / 2;
static constexpr const double dnmt_GAUSS_XUPPER = 8.572;
static constexpr const double dnmt_GAUSS_XLOWER = -37.519;
static constexpr const double dnmt_GAUSS_SCALE = 16.0;

[[nodiscard]] static double
get_del(const double x, const double rational) {
  static constexpr auto coeff = -0.5;
  const double xsq = std::floor(x * dnmt_GAUSS_SCALE) / dnmt_GAUSS_SCALE;
  const double del = (x - xsq) * (x + xsq);
  return std::exp(coeff * (xsq * xsq + del)) * rational;
}

// Normal cdf for fabs(x) < 0.66291
[[nodiscard]] static double
gauss_small(const double x) {
  // clang-format off
  static constexpr auto a_last = 4;
  static constexpr auto a = std::array{
    2.2352520354606839287,
    161.02823106855587881,
    1067.6894854603709582,
    18154.981253343561249,
    0.065682337918207449113,
  };
  static constexpr auto b_last = 3;
  static constexpr auto b = std::array{
    47.20258190468824187,
    976.09855173777669322,
    10260.932208618978205,
    45507.789335026729956,
  };
  // clang-format on

  const double xsq = x * x;
  double xnum = a[a_last] * xsq;
  double xden = xsq;
  for (std::uint32_t i = 0; i < b_last; i++) {
    xnum = (xnum + a[i]) * xsq;  // NOLINT(*-constant-array-index)
    xden = (xden + b[i]) * xsq;  // NOLINT(*-constant-array-index)
  }

  return x * (xnum + a[b_last]) / (xden + b[b_last]);
}

// Normal cdf for 0.66291 < fabs(x) < sqrt(32).
[[nodiscard]] static double
gauss_medium(const double x) {
  static constexpr auto c_last = 8;
  static constexpr auto d_last = 7;
  // clang-format off
  static constexpr auto c = std::array{
    0.39894151208813466764,
    8.8831497943883759412,
    93.506656132177855979,
    597.27027639480026226,
    2494.5375852903726711,
    6848.1904505362823326,
    11602.651437647350124,
    9842.7148383839780218,
    1.0765576773720192317e-8,
  };
  static constexpr auto d = std::array{
    22.266688044328115691,
    235.38790178262499861,
    1519.377599407554805,
    6485.558298266760755,
    18615.571640885098091,
    34900.952721145977266,
    38912.003286093271411,
    19685.429676859990727,
  };
  // clang-format on
  const double absx = std::abs(x);

  double xnum = c[c_last] * absx;
  double xden = absx;

  for (unsigned int i = 0; i < d_last; i++) {
    xnum = (xnum + c[i]) * absx;  // NOLINT(*-constant-array-index)
    xden = (xden + d[i]) * absx;  // NOLINT(*-constant-array-index)
  }

  const double temp = (xnum + c[d_last]) / (xden + d[d_last]);

  return get_del(x, temp);
}

// Normal cdf for
// { sqrt(32) < x < GAUSS_XUPPER } union { GAUSS_XLOWER < x < -sqrt(32) }
[[nodiscard]] static double
gauss_large(const double x) {
  static constexpr auto p_last = 5;
  static constexpr auto q_last = 4;
  // clang-format off
  static constexpr auto p = std::array{
    0.21589853405795699,
    0.1274011611602473639,
    0.022235277870649807,
    0.001421619193227893466,
    2.9112874951168792e-5,
    0.02307344176494017303,
  };
  static constexpr auto q = std::array{
    1.28426009614491121,
    0.468238212480865118,
    0.0659881378689285515,
    0.00378239633202758244,
    7.29751555083966205e-5,
  };
  // clang-format on

  const double absx = std::abs(x);
  const double xsq = 1.0 / (x * x);
  double xnum = p[p_last] * xsq;
  double xden = xsq;

  for (int i = 0; i < q_last; i++) {
    xnum = (xnum + p[i]) * xsq;  // NOLINT(*-constant-array-index)
    xden = (xden + q[i]) * xsq;  // NOLINT(*-constant-array-index)
  }

  double temp = xsq * (xnum + p[q_last]) / (xden + q[q_last]);
  temp = (dnmt_M_1_SQRT2PI - temp) / absx;

  return get_del(x, temp);
}

[[nodiscard]] double
dnmt_gsl_cdf_ugaussian_P(const double x) {
  static constexpr auto absx_cutoff = 0.66291;
  const double absx = std::abs(x);
  if (absx < dnmt_GAUSS_EPSILON)
    return 0.5;  // NOLINT(*-avoid-magic-numbers)
  if (absx < absx_cutoff)
    return 0.5 + gauss_small(x);  // NOLINT(*-avoid-magic-numbers)
  if (absx < dnmt_SQRT32) {
    const auto result = gauss_medium(x);
    return x > 0.0 ? 1.0 - result : result;
  }
  if (x > dnmt_GAUSS_XUPPER)
    return 1.0;
  if (x < dnmt_GAUSS_XLOWER)
    return 0.0;
  const auto result = gauss_large(x);
  return x > 0.0 ? 1.0 - result : result;
}

[[nodiscard]] double
dnmt_gsl_cdf_ugaussian_Q(const double x) {
  static constexpr auto absx_cutoff = 0.66291;
  const double absx = std::abs(x);
  if (absx < dnmt_GAUSS_EPSILON)
    return 0.5;  // NOLINT(*-avoid-magic-numbers)
  if (absx < absx_cutoff) {
    const auto result = gauss_small(x);
    // NOLINTNEXTLINE(*-avoid-magic-numbers)
    return x < 0.0 ? std::abs(result) + 0.5 : 0.5 - result;
  }
  if (absx < dnmt_SQRT32) {
    const auto result = gauss_medium(x);
    return x < 0.0 ? 1.0 - result : result;
  }
  if (x > -dnmt_GAUSS_XLOWER)
    return 0.0;
  if (x < -dnmt_GAUSS_XUPPER)
    return 1.0;
  const auto result = gauss_large(x);
  return x < 0.0 ? 1.0 - result : result;
}

[[nodiscard]] double
dnmt_gsl_cdf_gaussian_P(const double x, const double sigma) {
  return dnmt_gsl_cdf_ugaussian_P(x / sigma);
}

[[nodiscard]] double
dnmt_gsl_cdf_gaussian_Q(const double x, const double sigma) {
  return dnmt_gsl_cdf_ugaussian_Q(x / sigma);
}
