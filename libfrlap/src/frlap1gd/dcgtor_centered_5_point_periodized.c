/*   libfrlap
 *
 *   src/frlap1gd/dcgtor_centered_5_point_periodized.c
 *     Generator of periodized 5-point centered differences
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
#include <stdio.h>  // fprintf()
#include <stdlib.h> // abort()

int
frlap1gd_dcgtor_centered_5_point_periodized (double alpha,
                                             double h,
                                             size_t n,
                                             double MU[const static n])

{
  (void) fprintf(stderr,
      "frlap1qd_dcgtor_centered_5_point_periodized() not implemented\n");
  abort();
  return 0;
}
