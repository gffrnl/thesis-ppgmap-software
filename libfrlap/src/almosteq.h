/*   libfrlap
 *
 *   src/almosteq.h
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

#pragma once

#include <stdbool.h>

/*
 * REMARK:
 *   - maybe put the functions almosteq in a library
 */
bool almosteqf (     float  a,      float  b);
bool almosteq  (     double a,      double b);
bool almosteql (long double a, long double b);
