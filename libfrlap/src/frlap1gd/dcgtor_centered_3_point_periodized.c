/*   libfrlap
 *
 *   src/frlap1gd/dcgtor_centered_3_point_periodized.c
 *     Generator of periodized 3-point centered differences
 *     coefficients.
 *
 *   Copyright (C) 2023  Guilherme F. Fornel <gffrnl@gmail.com>
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
#include <frlap/frlap.h>
#include <frlap/frlap1gd/dcgtors.h>
#include <math.h>

/*
 * TODO:
 *   - maybe put the beta function in a library of
 *     special math functions
 */

/**
 *  Math special function beta(a, b);
 *    - here a, b are reals, but in fact can be complex
 *    - definition:
 *        beta(a, b) =
 *          (gamma(a) * gamma(b)) / gamma(a + b) 
 */
static double beta (double a, double b)
{
 /*
  * REMARK: using log-gamma function and then exponentiate
  *         is better to mitigate round-off/approximation
  *         errors
  */
  return exp(lgamma(a) + lgamma(b) - lgamma(a + b));
}

void
frlap1qd_dcgtor_centered_3_point_periodized (double frac_expon,
                                             double grid_step,
                                             size_t n,
                                             double mu[const static n])
{
  double c1 = 0.0, c2 = 0.0;

  c1 = 1.0 / (pow(grid_step, - frac_expon) * M_PI);
  c2 = c1 * sin(0.5 * frac_expon * M_PI);
  
  mu[0] = - c1 * exp2(frac_expon) * beta(0.5 + 0.5 * frac_expon, 0.5);

  for (size_t k = 1; k < n; k++)
    mu[k] = c2 * beta((double) k - 0.5 * frac_expon, 1.0 + frac_expon);
}
