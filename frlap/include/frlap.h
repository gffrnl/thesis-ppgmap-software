/* frlap.h */
#pragma once
#include <stddef.h>
#include <math.h>

double frlap_ncons (size_t n, double order);
double frlap_ncons1(double order);



/*
#include <stddef.h>
#include <assert.h>
#include <math.h>
#include "mu.h"
#include "convolve.h"

static
void
frlap1 (const double alpha, const double h,
        const size_t n,
        const double y[restrict const static n],
              double l[restrict const static n])
{
  double *mu = NULL;

  mu = (double *) malloc(n * sizeof(double));
  assert(mu != NULL);
  muper3(alpha, h, n, mu);
  for (size_t i = 0; i < n; i++)
    mu[i] *= -1.0;
  (void) symmconvolve(n, mu, y, l);
  free(mu);
}

static
void
tfrlap1 (const double alpha, const double h,
         const size_t n,
         const double y[restrict const static n],
               double l[restrict const static n])
{
  double *mu = NULL;

  mu = (double *) malloc(n * sizeof(double));
  assert(mu != NULL);
  muper3(alpha, h, n, mu);
  for (size_t i = 0; i < n; i++)
    mu[i] *= -1.0;
  (void) symmconvolve(n, mu, y, l);
  free(mu);
}
*/
