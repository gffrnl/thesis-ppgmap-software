/*   libfrlap
 *
 *   src/frlap1gd/dcgtor_gorenflo_mainardi.c
 *     Generator of Gorenflo & Mainardi differences coefficients.
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
#include <stdio.h>   // TODO: Maybe remove

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

void
frlap1qd_dcgtor_gorenflo_mainardi (double frac_expon,
                                   double grid_step,
                                   size_t n,
                                   double mu[const static n])
{
  /*
   * TODO:
   *   - seek for round-off errors, catastrophic
   *     cancelation or other numeric problems
   *   - maybe rewrite for optimization
   */

  double c1 = 0.0, c2 = 0.0;

  if (almost_equal(frac_expon, 1.0))
    c1 = 1.0 / (grid_step * M_PI);
  else
    {
      c1 = - pow(grid_step, - frac_expon)
        / ( 2.0 * cos(0.5 * grid_step * M_PI) );

      c2 = c1 * frac_expon * (frac_expon - 1)
        / tgamma(2.0 - frac_expon);
    }

  if (almost_equal(frac_expon, 1.0))
    {
      mu[0] = - 2.0 * c1;

      for (size_t k = 1; k < n; k++)
          mu[k] = c1 / (k * (k+1));
    }
  else if (frac_expon < 1.0)
    {
      mu[0] = 2.0 * c1;

      /*
       * REMARK: using log-gamma function and then
       *         exponentiate is better to mitigate
       *         round-off/approximation errors
       */
      for (size_t k = 1; k < n; k++)
        mu[k] = c2 * exp( lgamma(k-frac_expon) - lgamma(k+1) );
    }
  else if (frac_expon > 1.0)
    {
      mu[0] = - 2.0 * frac_expon * c1;

      mu[1] = c1 * (1.0 + 0.5 * frac_expon * (frac_expon - 1.0));

      /*
       * REMARK: using log-gamma function and then
       *         exponentiate is better to mitigate
       *         round-off/approximation errors
       */
      for (size_t k = 2; k < n; k++)
        mu[k] = c2 * exp( lgamma(k+1-frac_expon) - lgamma(k+2) );
    }
}
