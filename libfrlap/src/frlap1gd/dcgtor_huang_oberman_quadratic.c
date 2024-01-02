/*   libfrlap
 *
 *   src/frlap1gd/dcgtor_huang_oberman_quadratic.c
 *     Generator of Huang & Oberman quadratic differences
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
#include "../almosteq.h"

#ifndef M_1_SQRTPI
#define M_1_SQRTPI 0.564189583547756286948079451561
#endif

int
frlap1gd_dcgtor_huang_oberman_quadratic (double alpha,
                                         double h,
                                         size_t n,
                                         double MU[const static n])
{
  extern double d0G_alpha_ne_1 (double, size_t);
  extern double d1G_alpha_ne_1 (double, size_t);
  extern double d2G_alpha_ne_1 (double, size_t);
  extern double d0G_alpha_eq_1 (size_t);
  extern double d1G_alpha_eq_1 (size_t);
  extern double d2G_alpha_eq_1 (size_t);

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

    double ca = 0.0;
    double ch = 0.0;
    
    if ( almosteq(alpha, 1.0) )
      {
        ca = - FRLAP_C1ALPHA(1.0);
        ch = 1.0 / h;
      }
    else
      {
        ca = - FRLAP_C1ALPHA(alpha);
        ch = pow(h, -alpha);
      }
    
    /*
     * REMARK: using Log-Gamma function and then
     *         exponentiate seems better to mitigate
     *         round-off/approximation errors.
     *         Since lgamma() computes the log of
     *         the abs value of Gamma, we only can
     *         do this because the arguments of
     *         our Gammas are always positive.
     */
    MU[0] = ch * M_1_SQRTPI * exp2(alpha)
      * exp(   lgamma(0.5 * (alpha + 1.0))
             - lgamma(2.0 -  0.5 * alpha ) );

    ch *= ca;
    if ( almosteq(alpha, 1.0) )
      {
        MU[1] = ch * (   1.0
                       -       d2G_alpha_eq_1(1)
                       - 0.5 * d1G_alpha_eq_1(3)
                       - 1.5 * d1G_alpha_eq_1(1)
                       +       d0G_alpha_eq_1(3)
                       -       d0G_alpha_eq_1(1) );
        
        for (size_t k = 2; k < n; k+=2)
          MU[k] = ch * 2.0 * (   d1G_alpha_eq_1(k+1)
                               + d1G_alpha_eq_1(k-1)
                               - d0G_alpha_eq_1(k+1)
                               + d0G_alpha_eq_1(k-1) );

        for (size_t k = 3; k < n; k+=2)
          MU[k] = ch * ( - 0.5 * (   d1G_alpha_eq_1(k+2)
                                   + d1G_alpha_eq_1(k-2) )
                         - 3.0 * d1G_alpha_eq_1(k)
                         +       d0G_alpha_eq_1(k+2)
                         -       d0G_alpha_eq_1(k-2) );
      }
    else
      {
        MU[1] = ch * (   1.0 / (2.0 - alpha)
                       -       d2G_alpha_ne_1(alpha, 1)
                       - 0.5 * d1G_alpha_ne_1(alpha, 3)
                       - 1.5 * d1G_alpha_ne_1(alpha, 1)
                       +       d0G_alpha_ne_1(alpha, 3)
                       -       d0G_alpha_ne_1(alpha, 1) );
        
        for (size_t k = 2; k < n; k+=2)
          MU[k] = ch * 2.0 * (   d1G_alpha_ne_1(alpha, k+1)
                               + d1G_alpha_ne_1(alpha, k-1)
                               - d0G_alpha_ne_1(alpha, k+1)
                               + d0G_alpha_ne_1(alpha, k-1) );

        for (size_t k = 3; k < n; k+=2)
          MU[k] = ch * ( - 0.5 * (   d1G_alpha_ne_1(alpha, k+2)
                                   + d1G_alpha_ne_1(alpha, k-2) )
                         - 3.0 * d1G_alpha_ne_1(alpha, k)
                         +       d0G_alpha_ne_1(alpha, k+2)
                         -       d0G_alpha_ne_1(alpha, k-2) );
      }    
  }

  return 0;
}
