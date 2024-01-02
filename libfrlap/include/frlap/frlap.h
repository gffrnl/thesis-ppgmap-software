/*   libfrlap
 *
 *   frlap.h - commmon macros and function definitions
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

double frlap_normal_const(unsigned n, double alpha);
double frlap_normal_const_1(double alpha);

#define FRLAP_CNALPHA(N, ALPHA) frlap_normal_const  (N, ALPHA)
#define FRLAP_C1ALPHA(ALPHA)    frlap_normal_const_1(ALPHA)
