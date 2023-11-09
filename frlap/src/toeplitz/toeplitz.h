/* toeplitz.h */
#pragma once
#include <stddef.h>

void
symm_toeplitz_prod (const size_t n,
    const double A1[const restrict static n],
    const double x [const restrict static n],
          double b [const restrict static n]);

int
fast_symm_toeplitz_prod (const size_t n,
    const double A1[const restrict static n],
    const double x [const restrict static n],
          double b [const restrict static n]);


void
durbin (const size_t n,
    const double r[restrict const static n],
          double y[restrict const static n]);

int
levinson (const size_t n,
    const double A1[restrict const static n],
          double x [restrict const static n],
    const double b [restrict const static n]);
    
int
levinson_in_place (const size_t n,
          double A1[restrict const static n],
          double x [restrict const static n],
          double b [restrict const static n]);
