/*   libfrlap
 *
 *   src/almosteq.c
 *     Comparison functions for floating point equality.
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

#include "almosteq.h"
#include <float.h>
#include <math.h>

/*
 * REMARK:
 *   - maybe put the functions almosteq in a library
 *   - see  https://en.cppreference.com/w/cpp/types/numeric_limits/epsilon
 *     for a better implementation of the equality between floating point
 *     types
 */

bool
almosteqf (float a, float b)
{
  float p = fabsf(a + b);
  float m = fabsf(a - b);

  if ( p < FLT_EPSILON )
    return true;

  return ( m < (FLT_EPSILON * p) );
}

bool
almosteq (double a, double b)
{
  double p = fabs(a + b);
  double m = fabs(a - b);

  if ( p < DBL_EPSILON )
    return true;

  return ( m < (DBL_EPSILON * p) );
}

bool
almosteql (long double a, long double b)
{
  long double p = fabsl(a + b);
  long double m = fabsl(a - b);

  if ( p < LDBL_EPSILON )
    return true;

  return ( m < (LDBL_EPSILON * p) );
}
