/*   libfrlap
 *
 *   frlap.c - commmon functions
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
#include <math.h>

#ifndef M_SQRTPI
#define M_SQRTPI  1.77245385090551602729816748334
#endif

double
frlap_normal_const (unsigned n, double alpha)
{
  // TODO: validate arguments
  return alpha * exp2(alpha - 1.0)
    / pow(M_SQRTPI, n)
    * exp(  lgamma( 0.5 * ((double) n + alpha) )
          - lgamma( 0.5 * (2.0        - alpha) ) );
}

double
frlap_normal_const_1 (double alpha)
{
  // TODO: validate arguments
  return alpha * exp2(alpha - 1.0) / M_SQRTPI
    * exp(  lgamma( 0.5 * (1.0 + alpha) )
          - lgamma( 0.5 * (2.0 - alpha) ) );
}
