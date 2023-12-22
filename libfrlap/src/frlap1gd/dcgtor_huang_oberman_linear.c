/*   libfrlap
 *
 *   src/frlap1gd/dcgtor_huang_oberman_linear.c
 *     Generator of Huang & Oberman linear differences
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
#include <float.h>

#include <stdbool.h> // TODO: Maybe remove

/*
 * TODO:
 *   - maybe put the function almost equal in a library
 *   - see  https://en.cppreference.com/w/cpp/types/numeric_limits/epsilon
 *     for a better implementation of the equality between floating point
 *     types
 */

static bool almost_equal(double a, double b)
{
  return fabs(a - b) < DBL_EPSILON;
}

/* TODO: REMOVE
static double frlap_ncons(size_t n, double order)
{
  return order * exp2(order - 1)
    / pow(M_PI, 0.5 * n)
      * exp(  lgamma(0.5 * (order + n))
            - lgamma( 1.0 - 0.5*order )  );
}
*/

/* TODO: REMOVE
static double frlap_ncons1(double order)
{
  return frlap_ncons(1, order);
}
*/

void
frlap1qd_dcgtor_huang_oberman_linear (double frac_expon,
                                      double grid_step,
                                      size_t n,
                                      double mu[const static n])
{
  extern double d1G_alpha_ne_1 (double, size_t);
  extern double d2G_alpha_ne_1 (double, size_t);

  extern double d1G_alpha_eq_1 (size_t);
  extern double d2G_alpha_eq_1 (size_t);

  double c1 = 0.0, c2 = 0.0;

  if (almost_equal(frac_expon, 1.0))
    c1 = 1.0 / grid_step;
  else
    c1 = pow(grid_step, - frac_expon);

  c2 = FRLAP_C1ALPHA(frac_expon) * c1;
  
 /*
  * REMARK: using log-gamma function and then
  *         exponentiate is better to mitigate
  *         round-off/approximation errors
  */
  mu[0] = - c1 * exp2(frac_expon) / sqrt(M_PI)
    * exp(  lgamma(0.5 * (1.0 + frac_expon))
          - lgamma(2.0 - 0.5 * frac_expon  ) );

  if (almost_equal(frac_expon, 1.0))
    {
      mu[1] = c2 * ( 1.0 - d2G_alpha_eq_1(1)
                         + d1G_alpha_eq_1(2)
                         - d1G_alpha_eq_1(1) );

      for (size_t k = 2; k < n; k++)
        mu[k] =
          c2 * (         d1G_alpha_eq_1(k+1)
                 +       d1G_alpha_eq_1(k-1)
                 - 2.0 * d1G_alpha_eq_1(k)   );
    }
  else
    {
      mu[1] = c2 * (   1.0 / (2.0 - frac_expon)
                     - d2G_alpha_ne_1(frac_expon, 1)
                     + d1G_alpha_ne_1(frac_expon, 2)
                     - d1G_alpha_ne_1(frac_expon, 1) );

      for (size_t k = 2; k < n; k++)
        mu[k] =
          c2 * (         d1G_alpha_ne_1(frac_expon, k+1)
                 +       d1G_alpha_ne_1(frac_expon, k-1)
                 - 2.0 * d1G_alpha_ne_1(frac_expon, k)   );
    }
}
