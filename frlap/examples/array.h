/* array.h */
#pragma once
#include <stdlib.h>
#include <stdio.h>

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
matmul (const size_t m,
        const size_t n,
        const double a[const static m][n],
        const double x[const static n],
              double b[const static m])
{
  for (size_t i = 0; i < m ; i++)
    {
      b[i] = 0.0;
      for (size_t j = 0; j < n; j++)
        b[i] += a[i][j] * x[j];
    }
}
