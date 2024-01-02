/*   libfrlap
 *
 *   src/frlap1gd/dcgtor_gorenflo_mainardi.c
 *     Generator for Gorenflo & Mainardi differences coefficients.
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

int
frlap1gd_dcgtor_gorenflo_mainardi (double alpha,
                                   double h,
                                   size_t n,
                                   double MU[const static n])
{
  /* TODO:
   *   - there is a way to mitigate the large errors
   *     when alpha ~~ 1 but alpha != 1 ???
   *   - when alpha == 1 it is better use 1.0/(k*(k+1))
   *     or take the logs and then exponentiate ???
   */
  
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
    if ( almosteq(alpha, 1.0) )
      {
        double c1 = 1.0 / (h * M_PI);
        
        MU[0] = c1 * 2.0;

        for (size_t k = 1; k < n; k++)
          MU[k] = - c1 / (k * (k+1));
          /* or it is better to use ??
          MU[k] = - c1 / exp(  log((double)  k   )
                             + log((double) (k+1)) );
          */
      }
    else
      {
        double ch = 1.0 / pow(h, alpha);
        double cc = 1.0 / cos(0.5 * alpha * M_PI);
        double ca = alpha;

        /*
         * REMARK: Using Log-Gamma function and then
         *         exponentiate seems to be better to
         *         mitigate round-off/approximation errors.
         *         We only can do this because the result
         *         of our Gammas are positive, since lgamma()
         *         computes the logarithm of Gamma's ABSOLUTE
         *         value.
         */

        if (alpha < 1.0)
          {
            MU[0] = cc * ch;
            
            ca *= 0.5 * (alpha - 1.0) * ch;
            for (size_t k = 1; k < n; k++)
              MU[k] = cc * ( ca * exp(   lgamma((double) k - alpha)
                                       - lgamma(k + 1)
                                       - lgamma(2.0 - alpha)  ) );
          }
        else /* 1 < alpha < 2 */
          {
            MU[0] = - cc * (ca * ch);

            ca *= 0.5 * (alpha - 1.0);
            MU[1] = cc * ( (0.5 + 0.5 * ca) * ch);

            ca *= (2.0 - alpha) * ch;
            for (size_t k = 2; k < n; k++)
              MU[k] = cc * ( ca * exp(   lgamma((double) (k+1) - alpha)
                                       - lgamma(k + 2)
                                       - lgamma(3.0 - alpha)  ) );
          }
      }
  }

  return 0;
}
