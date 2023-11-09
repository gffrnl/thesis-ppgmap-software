/* frlap/frlap1dq.h */
#pragma once
#include <stddef.h>
#include <stdlib.h>
#include <stdbool.h> // needed before C23
#include <frlap/frlap1dq/wgtor.h>

#include <assert.h> // TODO: Maybe delete this line
                    //       (if not to use assert)

#include <stdio.h>
//#include "../array.h"

#include <toeplitz/toeplitz.h>

enum frlap1dq_wtype {
  FRLAP1DQ_WTYPE_SPECTRAL,
  FRLAP1DQ_WTYPE_HUANG_OBERMAN_LINEAR,
  FRLAP1DQ_WTYPE_HUANG_OBERMAN_QUADRATIC,
  FRLAP1DQ_WTYPE_GORENFLO_MAINARDI,
  FRLAP1DQ_WTYPE_CENTERED_PERIODIC_3_POINT,
  FRLAP1DQ_WTYPE_CENTERED_PERIODIC_5_POINT,
};

struct frlap1dq {
  double order;
  double gstep;
  enum frlap1dq_wtype wtype;
  void (*wgtor)(double, double, size_t n, double[const static n]);
};

void
tfrlap1dq_g (struct frlap1dq prm,
             const size_t n,
             const double y[const static n],
                   double l[const static n])
{
  assert(prm.wgtor != NULL);
  double *w = NULL;

  // Generation of the weigths:
  w = (double *) malloc(n * sizeof(double));
  assert(w != NULL);
  prm.wgtor(prm.order, prm.gstep, n, w);

  fast_symm_toeplitz_prod(n, w, y, l);

  free(w);

  for (size_t i = 0; i < n; i++)
    l[i] *= -1.0;
}

void
tfrlap1dq_g__naive__ (struct frlap1dq prm,
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

  // puts("\n ALLOCATED SPACE FOR THE WEIGHTS");
    

  prm.wgtor(prm.order, prm.gstep, n, w);
  //puts("\nw:");
  //print_arr(n, w);

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

void
tfrlap1dq (struct frlap1dq prm,
           const size_t n,
           const double y[const restrict static n],
                 double l[const restrict static n])
{
  switch (prm.wtype)
    {
      /*
    case FRLAP1DQ_WTYPE_SPECTRAL:
      prm.wgtor = &frlap1dq_wgtor_spectral;
      break;
      */
    case FRLAP1DQ_WTYPE_HUANG_OBERMAN_LINEAR:
      prm.wgtor = &frlap1dq_wgtor_huang_oberman_linear;
      //puts("\n USING WEIGHTS HUANG-OBERMAN LINEAR");
      break;
    case FRLAP1DQ_WTYPE_HUANG_OBERMAN_QUADRATIC:
      prm.wgtor = &frlap1dq_wgtor_huang_oberman_quadratic;
      //puts("\n USING WEIGHTS HUANG-OBERMAN QUADRATIC");
      break;
      /*
    case FRLAP1DQ_WTYPE_GORENFLO_MAINARDI:
      prm.wgtor = &frlap1dq_wgtor_gorenflo_mainardi;
      break;
    case FRLAP1DQ_WTYPE_CENTERED_PERIODIC_3_POINT:
      prm.wgtor = &frlap1dq_wgtor_centered_periodic_3_point;
      break;
    case FRLAP1DQ_WTYPE_CENTERED_PERIODIC_5_POINT:
      prm.wgtor = &frlap1dq_wgtor_centered_periodic_5_point;
      break;
      */

    default:
      assert(("weight type not known", false));
      return;
    }
  
  tfrlap1dq_g(prm, n, y, l);
}

void
tfrlap1dq__naive__ (struct frlap1dq prm,
                    const size_t n,
                    const double y[const restrict static n],
                          double l[const restrict static n])
{
  switch (prm.wtype)
    {
      /*
    case FRLAP1DQ_WTYPE_SPECTRAL:
      prm.wgtor = &frlap1dq_wgtor_spectral;
      break;
      */
    case FRLAP1DQ_WTYPE_HUANG_OBERMAN_LINEAR:
      prm.wgtor = &frlap1dq_wgtor_huang_oberman_linear;
      //puts("\n USING WEIGHTS HUANG-OBERMAN LINEAR");
      break;
    case FRLAP1DQ_WTYPE_HUANG_OBERMAN_QUADRATIC:
      prm.wgtor = &frlap1dq_wgtor_huang_oberman_quadratic;
      //puts("\n USING WEIGHTS HUANG-OBERMAN QUADRATIC");
      break;
      /*
    case FRLAP1DQ_WTYPE_GORENFLO_MAINARDI:
      prm.wgtor = &frlap1dq_wgtor_gorenflo_mainardi;
      break;
    case FRLAP1DQ_WTYPE_CENTERED_PERIODIC_3_POINT:
      prm.wgtor = &frlap1dq_wgtor_centered_periodic_3_point;
      break;
    case FRLAP1DQ_WTYPE_CENTERED_PERIODIC_5_POINT:
      prm.wgtor = &frlap1dq_wgtor_centered_periodic_5_point;
      break;
      */

    default:
      assert(("weight type not known", false));
      return;
    }
  
  tfrlap1dq_g__naive__(prm, n, y, l);
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------
//--------------------------------------------------------------------
//--------------------------------------------------------------------


void
tfrlap1dq_0g (struct frlap1dq prm,
              const size_t n,
              const double y[const static n],
              const size_t ja,
              const size_t jb,
                    double l[const static jb-ja+1])
{
  double *mu = NULL;
  double *la = NULL;
  double *lb = NULL;
  double (*Ba)[ja] = NULL;
  double (*Bb)[(n-1) - (ptrdiff_t) jb] = NULL;
  size_t m = 0;
  
  // Check the ranges
  assert(ja >= 0);
  assert(jb > ja);
  assert(jb < n);

  
  // Check the weight generator
  assert(prm.wgtor != NULL);

  m = jb-ja+1;

  /*
  (void) printf("ja = %lld, jb = %lld, m = %lld\n\n",
                ja, jb, m);
  */
  
  // Generation of the weigths:
  mu = (double *) malloc(n * sizeof(double));
  assert(mu != NULL);
  prm.wgtor(prm.order, prm.gstep, n, mu);


  assert(m > 0);

  /*
  for(size_t i = 0; i < n; i++)
    printf("mu[%d] = %16.8e\n", i, mu[i]);
  puts("\n");
  */
  /*
  for(size_t i = ja; i < ja+m; i++)
    printf("y[%d] = %16.8e\n", i, y[i]);
  puts("\n");
  */
  
  // Computation of the interior term:
  fast_symm_toeplitz_prod(m, mu, y+ja, l);

  /*
  for(size_t i = 0; i < m; i++)
    printf("l[%d] = %16.8e\n", i, l[i]);
  puts("\n");
  */
  
  // Computation of the boundary term:
  { 
    //  1. Creation of the 2d array that represents
    //     the matrix Ba
    Ba = (double(*)[ja]) calloc(m, sizeof(double[ja]));
    assert(Ba != NULL);

    for (size_t i = 0; i < m; i++)
      for (ptrdiff_t j = 0; j < ja; j++)
        Ba[i][j] = mu[((ptrdiff_t) (ja+i))-j];

    /*
    (void) puts("Matrix Ba:");
    matprint(m, ja, Ba);
    puts("\n");
    */
    
    la = (double *) malloc(m * sizeof(double));

    /*
    for (size_t i = 0; i < m; i++)
      for (ptrdiff_t j = 0; j < ja; j++)
        la[i] = Ba[i][j] * y[j];
    */

    matmul(m, ja, Ba, y, la);

    /*
    for(size_t i = 0; i < m; i++)
      printf("la[%d] = %16.8e\n", i, la[i]);
    puts("\n");
    */
    
    free(Ba);

    for (size_t i = 0; i < m; i++)
      l[i] += la[i];

    free(la);


    //  2. Creation of the 2d array that represents
    //     the matrix Bb
    Bb = (double(*)[(n-1) - (ptrdiff_t) jb])
      calloc(m, sizeof(double[(n-1) - (ptrdiff_t) jb]));
    assert(Bb != NULL);

    /*
    (void) printf("(n-1) - (ptrdiff_t) jb = %lld\n\n",
                  (size_t) ((n-1) - (ptrdiff_t) jb));
    */
                  
    for (ptrdiff_t i = 0; i < m; i++)
      for (size_t j = 0; j < (n-1) - (ptrdiff_t) jb; j++)
        Bb[i][j] = mu[((ptrdiff_t) (m+j)) - i];

    /*
    (void) puts("Matrix Bb:");
    matprint(m, (n-1) - (ptrdiff_t) jb, Bb);
    puts("\n");
    */
    
    lb = (double *) malloc(m * sizeof(double));

    /*
    for (size_t i = 0; i < m; i++)
      for (size_t j = jb+1; j < n; j++)
        lb[i] = Bb[i][j - (ptrdiff_t) (jb+1)] * y[j];
    */
    matmul(m, (n-1) - (ptrdiff_t) jb, Bb, y+jb+1, lb);

    /*
    for(size_t i = 0; i < m; i++)
      printf("lb[%d] = %16.8e\n", i, lb[i]);
    puts("\n");
    */
    
    free(Bb);

    for (size_t i = 0; i < m; i++)
      l[i] += lb[i];

    free(lb);
  }
  
  free(mu);

  for (size_t i = 0; i < m; i++)
    l[i] *= -1.0;
}


void
tfrlap1dq_0 (struct frlap1dq prm,
             const size_t n,
             const double y[const restrict static n],
             const size_t ja,
             const size_t jb,
                   double l[const restrict static jb-ja+1])
{
  switch (prm.wtype)
    {
      /*
    case FRLAP1DQ_WTYPE_SPECTRAL:
      prm.wgtor = &frlap1dq_wgtor_spectral;
      break;
      */
    case FRLAP1DQ_WTYPE_HUANG_OBERMAN_LINEAR:
      prm.wgtor = &frlap1dq_wgtor_huang_oberman_linear;
      //puts("\n USING WEIGHTS HUANG-OBERMAN LINEAR");
      break;
    case FRLAP1DQ_WTYPE_HUANG_OBERMAN_QUADRATIC:
      prm.wgtor = &frlap1dq_wgtor_huang_oberman_quadratic;
      //puts("\n USING WEIGHTS HUANG-OBERMAN QUADRATIC");
      break;
      /*
    case FRLAP1DQ_WTYPE_GORENFLO_MAINARDI:
      prm.wgtor = &frlap1dq_wgtor_gorenflo_mainardi;
      break;
    case FRLAP1DQ_WTYPE_CENTERED_PERIODIC_3_POINT:
      prm.wgtor = &frlap1dq_wgtor_centered_periodic_3_point;
      break;
    case FRLAP1DQ_WTYPE_CENTERED_PERIODIC_5_POINT:
      prm.wgtor = &frlap1dq_wgtor_centered_periodic_5_point;
      break;
      */

    default:
      assert(("weight type not known", false));
      return;
    }
  
  tfrlap1dq_0g(prm, n, y, ja, jb, l);
}
