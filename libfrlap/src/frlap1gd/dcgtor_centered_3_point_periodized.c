/*   libfrlap
 *
 *   src/frlap1gd/dcgtor_centered_3_point_periodized.c
 *     Generator for 3-point centered periodized differences
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

#include <frlap/frlap1gd/dcgtors.h>
#include <math.h>
#include "../almosteq.h"

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
 // !!! a and b must be positive !!!
 /*
  * REMARK: using log-gamma function and then exponentiate
  *         is better to mitigate round-off/approximation
  *         errors
  */
  return exp(lgamma(a) + lgamma(b) - lgamma(a + b));
}

int
frlap1gd_dcgtor_centered_3_point_periodized (double alpha,
                                             double h,
                                             size_t n,
                                             double MU[const static n])
{
  { /* check arguments */
    int ret;
    ret = dcgtors_validate_args(alpha, h, n);
    if ( ret != 0 )
      return ret;
  }

  { /* treat the boundary cases alpha = 0 and alpha = 2 */
    if ( almosteq(alpha, 0.0) )
      {
        MU[0] = 1.0;
        for (size_t k = 1; k < n; k++)
          MU[k] = 0.0;
        return 0;
      }

    if ( almosteq(alpha, 2.0) )
      {
        MU[0] =  2.0 / (h * h);
        MU[1] = -1.0 / (h * h);
        for (size_t k = 2; k < n; k++)
          MU[k] = 0.0;
        return 0;
      }
  }
  
  { /* {0 < alpha < 2} */
    double ch = 1.0 / ( M_PI * pow(h, alpha) );
    double cs = - ch * sin(0.5 * alpha * M_PI);

    MU[0] = ch * ( exp2(alpha) * beta(0.5 + 0.5 * alpha, 0.5) );

    for (size_t k = 1; k < n; k++)
      MU[k] = cs * beta((double) k - 0.5 * alpha, 1.0 + alpha);
  }

  return 0;
}
