/*   libfrlap
 *
 *   include/frlap/frlap1gd/dcgtors_validate_args.h
 *     Argument validator for differences coefficients generators.
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
#pragma once
#include <stddef.h>
//#include <frlap/frlap1gd/dcgtors_errno.h>

#define FRLAP_FRLAP1GD_DCGTORS_ERRNO_ARG_FRAC_EXPON 1
#define FRLAP_FRLAP1GD_DCGTORS_ERRNO_ARG_GRID_STEP  2
#define FRLAP_FRLAP1GD_DCGTORS_ERRNO_ARG_COEFF_SIZE 3

int
dcgtors_validate_args (double frac_expon,
                       double grid_step,
                       size_t n);
