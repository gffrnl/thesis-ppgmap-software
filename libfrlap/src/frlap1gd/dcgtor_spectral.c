/*   libfrlap
 *
 *   src/frlap1gd/dcgtor_spectral.c
 *     Generator of spectral differences coefficients.
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
#include <ooura/intde/intde1.c>
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

static int freq;
static double expon;

static double f(double x)
{
  extern int freq;
  extern double expon;

  return pow(x, expon) * cos(freq * x);
}

void
frlap1qd_dcgtor_spectral (double frac_expon,
                          double grid_step,
                          size_t n,
                          double mu[const static n])
{
  if (almost_equal(frac_expon, 1.0))
    {
      mu[0] = -1.0;
      for (size_t k = 1; k < n; k+=2)
          mu[k] = 2.0 / (M_PI * k*k * grid_step);
      for (size_t k = 2; k < n; k+=2)
          mu[k] = 0.0;
    }
  else if (almost_equal(frac_expon, 2.0))
    {
      mu[0] = - (M_PI / grid_step) * (M_PI / grid_step) / 3.0;
      for (size_t k = 1; k < n; k+=2)
          mu[k] =   2.0 / ( (k*grid_step) * (k*grid_step) );
      for (size_t k = 2; k < n; k+=2)
          mu[k] = - 2.0 / ( (k*grid_step) * (k*grid_step) );
    }
  else
    {
      extern int freq;
      extern double expon;
      double err;

      expon = frac_expon;
      
      mu[0] = - pow(M_PI / grid_step, frac_expon) / (frac_expon + 1.0);

      grid_step = pow(grid_step, frac_expon);

      for (size_t k = 1; k < n; k++)
        {
          freq = k;
          intde(f, 0.0, M_PI, 1.1e-14, &mu[k], &err);
          mu[k] *= -1.0 / (M_PI * grid_step);
        }
    }
}
