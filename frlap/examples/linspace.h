/* linspace.h */
#pragma once
#include <stddef.h>
#include <stdbool.h>

enum linspace {
  EXCLUDE_ENDPOINT = 0,
  INCLUDE_ENDPOINT = 1
};

static
double
linspace_e (const double  a, const double b,
            enum linspace opt,
            const size_t  n,
                  double  x[const static n])
{
  double h = (b - a) /
      ( (opt == INCLUDE_ENDPOINT) ? (n - 1) : n );
  for (size_t i = 0; i < n; i++)
      x[i] = a + i * h;
  if (opt == INCLUDE_ENDPOINT)
      x[n-1] = b;
  return h;
}

static
double
linspace (const double a, const double  b,
          const size_t n,
                double x[const static n])
{
  return linspace_e(a, b, INCLUDE_ENDPOINT, n, x);
}
