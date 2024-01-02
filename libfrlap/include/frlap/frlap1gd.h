/*   libfrlap
 *
 *   frlap1gd.h - generalized differences for frac. Laplacian in dim. 1
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
//#include <stdlib.h>
//#include <stdbool.h> // needed before C23
//#include "wgtor/wgtor.h"

//#include <assert.h> // TODO: Maybe delete this line
                    //       (if not to use assert)

//#include <stdio.h>
//#include "../array.h"

enum frlap1gd_dctype {
  FRLAP1GD_DCTYPE_UNKNOWN,
  FRLAP1GD_DCTYPE_SPECTRAL,
  FRLAP1GD_DCTYPE_HUANG_OBERMAN_LINEAR,
  FRLAP1GD_DCTYPE_HUANG_OBERMAN_QUADRATIC,
  FRLAP1GD_DCTYPE_GORENFLO_MAINARDI,
  FRLAP1GD_DCTYPE_CENTERED_PERIODIC_3_POINT,
  FRLAP1GD_DCTYPE_CENTERED_PERIODIC_5_POINT,
};

struct frlap1gd {
  double frexpo;
  double step;
  enum frlap1gd_dctype dctype;
  void (*dcgtor)(double, double, size_t n, double[const static n]);
};



/*
void
frlap1gd_set_params(*struct frlap1gd,
                    alpha, 0.1,
                    FRLAP1GD_DCTYPE_HUANG_OBERMAN_QUADRATIC);
*/

/*

void
frlap1gd_gen_dc (struct frlap1gd prm,
                 size_t          n,
                 double          mu[const static n]);


size_t
frlap1gd_gen_dc (struct frlap1gd prm,
                 double          mu[const static n]);

*/

#define FRLAP1GD0_DIFF_COEFF_SIZE(n, ja, jb) (jb+1>n-ja)?jb+1:n-ja
//#define FRLAP1GD_FRAC_LAP_SIZE(n, ja, jb) jb-ja+1
#define FRLAP1GD0_FRLAP_SIZE(n, ja, jb) jb-ja+1

/*
void
tfrlap1gd0(size_t n,
           const double Y [static const n],
           size_t ja, size_t jb,
           size_t nc, // must be FRLAP1GD_DIFF_COEFF_SIZE(n, ja, jb)
           const double mu[static const nc],
           size_t n0, // must be FRLAP1GD_RETURN_SIZE(n, ja, jb)
           double FLY[static const n0]);

/*
void
tfrlap1qd_g (struct frlap1qd prm,
             const size_t n,
             const double y[const restrict static n],
                   double l[const restrict static n])
{
  assert(("function tfrlap1qd_g not implemented", false));
}
*/

/*
void
tfrlap1qd_g__naive__ (struct frlap1qd prm,
                      const size_t n,
                      const double y[const restrict static n],
                            double l[const restrict static n])
{
  double *w = NULL;
  
  assert(prm.wgtor != NULL);

  // TODO: generate the weights!

  // Generation of the weigths:
  w = (double *) malloc(n * sizeof(double));
  assert(w != NULL);

  puts("\n ALLOCATED SPACE FOR THE WEIGHTS");
    

  prm.wgtor(prm.order, prm.gstep, n, w);
  puts("\nw:");
  print_arr(n, w);

  for (size_t j = 0; j < n; j++)
    {
      l[j] = 0;
      for (size_t k = 0; k < n; k++)
        {
          l[j] -= w[abs(j-k)] * y[k];
        }
    }
    
  free(w);
}
*/

/*
void
tfrlap1qd (struct frlap1qd prm,
           const size_t n,
           const double y[const restrict static n],
                 double l[const restrict static n])
{
  assert(("function tfrlap1qd not implemented", false));
}
*/

/*
void
tfrlap1qd__naive__ (struct frlap1qd prm,
                    const size_t n,
                    const double y[const restrict static n],
                          double l[const restrict static n])
{
  //assert(("function tfrlap1qd__naive__ not implemented", false));

  switch (prm.wtype)
    {
*/
      /*
    case FRLAP1QD_WTYPE_SPECTRAL:
      prm.wgtor = &frlap1qd_wgtor_spectral;
      break;
      */
/*
    case FRLAP1QD_WTYPE_HUANG_OBERMAN_LINEAR:
      prm.wgtor = &frlap1qd_wgtor_huang_oberman_linear;
      puts("\n USING WEIGHTS HUANG-OBERMAN LINEAR");
      break;
    case FRLAP1QD_WTYPE_HUANG_OBERMAN_QUADRATIC:
      prm.wgtor = &frlap1qd_wgtor_huang_oberman_quadratic;
      puts("\n USING WEIGHTS HUANG-OBERMAN QUADRATIC");
      break;
*/
      /*
    case FRLAP1QD_WTYPE_GORENFLO_MAINARDI:
      prm.wgtor = &frlap1qd_wgtor_gorenflo_mainardi;
      break;
    case FRLAP1QD_WTYPE_CENTERED_PERIODIC_3_POINT:
      prm.wgtor = &frlap1qd_wgtor_centered_periodic_3_point;
      break;
    case FRLAP1QD_WTYPE_CENTERED_PERIODIC_5_POINT:
      prm.wgtor = &frlap1qd_wgtor_centered_periodic_5_point;
      break;
      */
/*
    default:
      assert(("weight type not known", false));
      return;
    }
  
  tfrlap1qd_g__naive__(prm, n, y, l);
}
*/
