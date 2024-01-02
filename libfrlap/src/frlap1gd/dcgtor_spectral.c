/*   libfrlap
 *
 *   src/frlap1gd/dcgtor_spectral.c
 *     Generator for spectral differences coefficients.
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

#include <ooura/intde/intde1.c>

static int freq;
static double expon;

static double f(double x)
{
  extern int freq;
  extern double expon;

  return pow(x, expon) * cos(freq * x);
}

int
frlap1gd_dcgtor_spectral (double alpha,
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
        double c2 = 2.0 / (h * h);
        
        MU[0] =  (M_PI * M_PI) / (3.0 * h * h);
        for (size_t k = 1; k < n; k+=2)
          MU[k] = - c2 / (k * k);
        for (size_t k = 2; k < n; k+=2)
          MU[k] =   c2 / (k * k);
        return 0;
      }
  }

  { /* the case alpha = 1 */
    if ( almosteq(alpha, 1.0) )
      {
        double c1 = 1.0 / (M_PI * h);
          
        MU[0] =  M_PI / (2 * h);
        for (size_t k = 1; k < n; k+=2)
          MU[k] = - c1 * (2.0 / (k * k));
        for (size_t k = 2; k < n; k+=2)
          MU[k] = 0.0;
        return 0;
      }
  }

  { /* {0 < alpha < 1} U {1 < alpha < 2} */
    extern int freq;
    extern double expon;
    double err;
    double ch = 1.0 / (M_PI * pow(h, alpha)); 
    
    MU[0] = pow(M_PI/h, alpha) / (alpha + 1.0);

    expon = alpha;
    for (size_t k = 1; k < n; k++)
      {
        freq = k;
        intde(f, 0.0, M_PI, 1.1e-14, &MU[k], &err);
        MU[k] *= ch;
      }
  }

  return 0;
}
