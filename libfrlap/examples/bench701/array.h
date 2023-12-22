/* array.h */
#pragma once
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <fltcmp.h>

static
size_t
seekval_dig(double val, size_t n, const double x[const static n], unsigned dig)
{
  size_t idx = n;
  for (size_t i = 0; i < n; i++)
    if (dblcmp_dig(x[i], val, dig))
      {
        idx = i;
        break;
      }
  return idx;
}

static
void
print_arr (const size_t n,
           const double x[const static n])
{
  (void) printf("[");
  for (size_t i = 0; i < n; i++)
    (void) printf("%.5e, ", x[i]);
  (void) puts("\b\b]");
}

static
void
elemental (const size_t   n,
           const double   x[restrict const static n],
                 double (*f) (double),
                 double   y[restrict const static n])
{
  for (size_t i = 0; i < n; i++)
      y[i] = f(x[i]);
}

static
void
matprint (const size_t m,
          const size_t n,
          const double a[const static m][n])
{
  for (size_t i = 0; i < m ; i++)
    {
      (void) printf("[");
      for (size_t j = 0; j < n; j++)
        {
          if (j == n-1)
            (void) printf("%.5e]\n", a[i][j]);
          else
            (void) printf("%.5e, ", a[i][j]);
        }
    }
}


static
void
matvecmul_sum (const size_t m,
               const size_t n,
               const double A[const static m][n],
               const double x[const static n],
                     double b[const static m])
{
  for (size_t i = 0; i < m ; i++)
    for (size_t j = 0; j < n; j++)
      b[i] += A[i][j] * x[j];
}

static
void
matvecmul (const size_t m,
           const size_t n,
           const double A[const static m][n],
           const double x[const static n],
                 double b[const static m])
{
  for (size_t i = 0; i < m ; i++)
    b[i] = 0.0;
  matvecmul_sum(m, n, A, x, b);
}
