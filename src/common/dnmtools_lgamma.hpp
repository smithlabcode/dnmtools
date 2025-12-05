/* Copyright (C) 2025 Andrew D. Smith
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

#ifndef DNMTOOLS_LGAMMA_HPP
#define DNMTOOLS_LGAMMA_HPP

/* ADS: the constants below, and some of the code, was taken from e_lgamma_r.c
   in Google's bionic source and lgamma_r.c in the MUSL source. One or both
   were modified from the original, maybe for style guidlines. Both files had
   this header.
 * ====================================================
 * Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
 *
 * Developed at SunSoft, a Sun Microsystems, Inc. business.
 * Permission to use, copy, modify, and distribute this
 * software is freely granted, provided that this notice
 * is preserved.
 * ====================================================
*/

#include <cassert>
#include <cmath>
#include <cstdint>

namespace dnmtools_lgamma {

// ADS: 'tc' is the unique minimum of the gamma function on R+
static constexpr std::double_t tc = +1.46163214496836224576e+00;
static constexpr std::double_t tf = -1.21486290535849611461e-01;
/* tt = -(tail of tf) */
static constexpr std::double_t tt = -3.63867699703950536541e-18;

// clang-format off

static constexpr std::double_t A0[6] = {
  +7.72156649015328655494e-02,
  +6.73523010531292681824e-02,
  +7.38555086081402883957e-03,
  +1.19270763183362067845e-03,
  +2.20862790713908385557e-04,
  +2.52144565451257326939e-05,
};

static constexpr std::double_t A1[6] = {
  +3.22467033424113591611e-01,
  +2.05808084325167332806e-02,
  +2.89051383673415629091e-03,
  +5.10069792153511336608e-04,
  +1.08011567247583939954e-04,
  +4.48640949618915160150e-05,
};

static constexpr std::double_t T0[5] = {
  +4.83836122723810047042e-01,
  -3.27885410759859649565e-02,
  +6.10053870246291332635e-03,
  -1.40346469989232843813e-03,
  +3.15632070903625950361e-04,
};

static constexpr std::double_t T1[5] = {
  -1.47587722994593911752e-01,
  +1.79706750811820387126e-02,
  -3.68452016781138256760e-03,
  +8.81081882437654011382e-04,
  -3.12754168375120860518e-04,
};

static constexpr std::double_t T2[5] = {
  +6.46249402391333854778e-02,
  -1.03142241298341437450e-02,
  +2.25964780900612472250e-03,
  -5.38595305356740546715e-04,
  +3.35529192635519073543e-04,
};

static constexpr std::double_t U[6] = {
  -7.72156649015328655494e-02,
  +6.32827064025093366517e-01,
  +1.45492250137234768737e+00,
  +9.77717527963372745603e-01,
  +2.28963728064692451092e-01,
  +1.33810918536787660377e-02,
};

static constexpr std::double_t V[5] = {
  +2.45597793713041134822e+00,
  +2.12848976379893395361e+00,
  +7.69285150456672783825e-01,
  +1.04222645593369134254e-01,
  +3.21709242282423911810e-03,
};

static constexpr std::double_t S[7] = {
  -7.72156649015328655494e-02,
  +2.14982415960608852501e-01,
  +3.25778796408930981787e-01,
  +1.46350472652464452805e-01,
  +2.66422703033638609560e-02,
  +1.84028451407337715652e-03,
  +3.19475326584100867617e-05,
};

static constexpr std::double_t R[6] = {
  +1.39200533467621045958e+00,
  +7.21935547567138069525e-01,
  +1.71933865632803078993e-01,
  +1.86459191715652901344e-02,
  +7.77942496381893596434e-04,
  +7.32668430744625636189e-06,
};

static constexpr std::double_t W[7] = {
  +4.18938533204672725052e-01,
  +8.33333333333329678849e-02,
  -2.77777777728775536470e-03,
  +7.93650558643019558500e-04,
  -5.95187557450339963135e-04,
  +8.36339918996282139126e-04,
  -1.63092934096575273989e-03,
};

// clang-format on

[[nodiscard]] inline auto
noneg_lgamma(double x) -> double {
  // ADS: this function does not work for very small or very large values.
  assert(std::isfinite(x) && x > 2.22045e-16);  // > 2^-52

  union {
    double f;  // unused
    std::uint64_t as_bits;
  } u = {x};  // safe inspection of binary layout

  // get magnitude
  const std::uint32_t ix = u.as_bits >> 32 & 0x7fffffff;

  std::double_t r{};
  // x exactly equal to 1 or 2
  if ((ix == 0x3ff00000 || ix == 0x40000000) &&
      static_cast<int>(u.as_bits) == 0)
    return 0;

  // ADS: for x < 2.22045e-16, just return the value of lgamma at
  // that 2.22045e-16

  // if (ix < 0x3CB00000)
  //   return 36.0437;

  // ADS: for x > 2^58 the code below applies. Maxing it out at 1.2994e+19,
  // because that's large enough and if we see this our estimation is going in
  // the wrong direction anyway.

  // if (ix > 0x43900000)
  //   return 1.12994e+19;  // proper way: return x * (std::log(x) - 1.0);

  // x < 2.0
  if (ix < 0x40000000) {
    int i{};
    std::double_t y{};
    if (ix <= 0x3feccccc) {  // lgamma(x) = lgamma(x+1)-log(x)
      r = -std::log(x);
      if (ix >= 0x3FE76944) {
        y = 1.0 - x;
        i = 0;
      }
      else if (ix >= 0x3FCDA661) {
        y = x - (tc - 1.0);
        i = 1;
      }
      else {
        y = x;
        i = 2;
      }
    }
    else {
      r = 0.0;
      if (ix >= 0x3FFBB4C3) {  // [1.7316,2]
        y = 2.0 - x;
        i = 0;
      }
      else if (ix >= 0x3FF3B4C4) {  // [1.23,1.73]
        y = x - tc;
        i = 1;
      }
      else {
        y = x - 1.0;
        i = 2;
      }
    }
    std::double_t p{};
    std::double_t w{}, z{};
    std::double_t p1{}, p2{}, p3{};
    switch (i) {
    case 0:
      z = y * y;
      p1 = A0[0] +
           z * (A0[1] + z * (A0[2] + z * (A0[3] + z * (A0[4] + z * A0[5]))));
      p2 =
        z * (A1[0] +
             z * (A1[1] + z * (A1[2] + z * (A1[3] + z * (A1[4] + z * A1[5])))));
      p = y * p1 + p2;
      r += -0.5 * y + p;
      break;
    case 1:
      z = y * y;
      w = z * y;
      p1 = T0[0] + w * (T0[1] + w * (T0[2] + w * (T0[3] + w * T0[4])));
      p2 = T1[0] + w * (T1[1] + w * (T1[2] + w * (T1[3] + w * T1[4])));
      p3 = T2[0] + w * (T2[1] + w * (T2[2] + w * (T2[3] + w * T2[4])));
      p = z * p1 - (tt - w * (p2 + y * p3));
      r += tf + p;
      break;
    case 2:
      p1 = y * (U[0] +
                y * (U[1] + y * (U[2] + y * (U[3] + y * (U[4] + y * U[5])))));
      p2 = 1.0 + y * (V[0] + y * (V[1] + y * (V[2] + y * (V[3] + y * V[4]))));
      r += -0.5 * y + p1 / p2;
    }
    return r;
  }

  // x < 8.0
  if (ix < 0x40200000) {
    const int i = static_cast<int>(x);
    const std::double_t y = x - static_cast<double>(i);
    const std::double_t p =
      y * (S[0] +
           y * (S[1] +
                y * (S[2] + y * (S[3] + y * (S[4] + y * (S[5] + y * S[6]))))));
    const std::double_t q =
      1.0 +
      y * (R[0] + y * (R[1] + y * (R[2] + y * (R[3] + y * (R[4] + y * R[5])))));
    r = 0.5 * y + p / q;    // RF P/Q
    std::double_t z = 1.0;  // lgamma(1+s) = log(s) + lgamma(s)
    // clang-format off
    switch (i) {
    case 7: z *= y + 6.0; // FALLTHRU
    case 6: z *= y + 5.0; // FALLTHRU
    case 5: z *= y + 4.0; // FALLTHRU
    case 4: z *= y + 3.0; // FALLTHRU
    case 3: z *= y + 2.0; // FALLTHRU
      r += std::log(z);
      break;
    }
    // clang-format on
    return r;
  }

  // 8.0 <= x (ix < 0x43900000)
  const std::double_t t = std::log(x);
  const std::double_t z = 1.0 / x;
  const std::double_t y = z * z;
  const std::double_t w =
    W[0] +
    z * (W[1] + y * (W[2] + y * (W[3] + y * (W[4] + y * (W[5] + y * W[6])))));
  return (x - 0.5) * (t - 1.0) + w;
}
};  // namespace dnmtools_lgamma

#endif  // DNMTOOLS_LGAMMA_HPP
