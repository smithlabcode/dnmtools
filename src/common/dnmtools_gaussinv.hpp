/* Code from GSl, see copyright below.
 */

/* cdf/gsl_cdf.h
 *
 * Copyright (C) 2002 Jason H. Stover.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  J. Stover */

double dnmt_gsl_cdf_ugaussian_Pinv(const double P);
double dnmt_gsl_cdf_ugaussian_Qinv(const double Q);

double dnmt_gsl_cdf_gaussian_P(const double x, const double sigma);
double dnmt_gsl_cdf_gaussian_Q(const double x, const double sigma);
