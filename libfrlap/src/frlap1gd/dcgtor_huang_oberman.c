/*   libfrlap
 *
 *   src/frlap1gd/dcgtor_huang_oberman.c
 *     Common functions used by Huang & Oberman differences
 *     coefficients generators.
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

#include <stddef.h>
#include <math.h>

double
d0G_alpha_ne_1 (double alpha, size_t k)
{
  return pow((double) k, 2.0 - alpha)
    / ( alpha * (alpha - 1.0) * (2.0 - alpha) );
}

double
d1G_alpha_ne_1 (double alpha, size_t k)
{
  return pow((double) k, 1.0 - alpha)
    / ( alpha * (alpha - 1.0) );
}

double
d2G_alpha_ne_1 (double alpha, size_t k)
{
  return - pow((double) k, -alpha) / alpha;
}

double
d0G_alpha_eq_1 (size_t k)
{
  return (double) k - k * log((double) k);
}

double
d1G_alpha_eq_1 (size_t k)
{
  return - log((double) k);
}

double
d2G_alpha_eq_1 (size_t k)
{
  return - 1.0 / (double) k;
}
