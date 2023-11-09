#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

double
static reverse_inner (n, x, y)
    const size_t n;
    const double x[restrict const static n];
    const double y[restrict const static n];
{
  size_t n1 = n - 1;
  double a = 0.0;
  for (size_t i = 0; i < n; ++i)
      a += x[n1-i] * y[i];
  return a;
}

int
levinson (n, A1, x, b)
    const size_t n;
    const double A1[restrict const static n];
          double x [restrict const static n];
    const double b [restrict const static n];
{
  size_t n1 = n - 1;
  double alpha, beta, mu, mult;
  double *r, *q, *y;

  /* Check data consistecy */
  if (fabs(A1[0]) < DBL_EPSILON)
      return 1; /* null diagonal! */

  /* Try allocate memory */
  r = (double *) malloc(n * sizeof(double));
  q = (double *) malloc(n * sizeof(double));
  y = (double *) malloc(n * sizeof(double));
  if (r == NULL || r == NULL || y == NULL)
    {
      /* cannot allocate memory:
         print to stderr and abort */
      (void) fprintf(stderr,
        "in function `levinson`: cannot "
        "allocate memory\n");
      abort();
    }

  for (size_t i = 0; i < n1; ++i)
      r[i] = A1[i+1] / A1[0];

  for (size_t i = 0; i < n; ++i)
      q[i] = b[i] / A1[0];

  y[0] = -r[0];
  x[0] =  q[0];
  alpha = y[0];
  mult = 1.0 - alpha*alpha;
  beta = 1.0;
  for (size_t k = 1; k < n; ++k)
    {
      size_t k1 = k-1;
      size_t k2 = (k+1)/2 - k%2;

      beta *= mult;
      mu = (q[k] - reverse_inner(k, r, x))
        / beta;
      for (size_t i = 0; i < k; ++i)
          x[i] += mu * y[k1-i];
      x[k] = mu;
      if (k < n1)
        { /* Durbin recursion */
          alpha =
            - (r[k] + reverse_inner(k, r, y))
            / beta;
          mult = 1.0 - alpha*alpha;
          for (size_t i = 0; i < k2; ++i)
            {
              y[i] += alpha * y[k1-i];
              y[k1-i] *= mult;
              y[k1-i] += alpha * y[i];
            }
          if (k%2 == 1)
              y[k2] *= (1.0 + alpha);
          y[k] = alpha;
        }        
    }
  free(q);
  free(y);
  free(r);
  return 0;
}

int
levinson_in_place (n, A1, x, b)
    const size_t n;
          double A1[restrict const static n];
          double x [restrict const static n];
          double b [restrict const static n];
{
  size_t n1 = n - 1;
  double alpha, beta, mu, mult;
  double A10 = A1[0];
  double *r = A1 + 1, *y;

  /* Check data consistecy */
  if (fabs(A10) < DBL_EPSILON)
      return 1; /* null diagonal! */

  /* Try allocate memory */
  y = (double *) malloc(n * sizeof(double));
  if (y == NULL)
    {
      /* cannot allocate memory:
         print to stderr and abort */
      (void) fprintf(stderr,
        "in function `levinson`: cannot "
        "allocate memory\n");
      abort();
    }

  for (size_t i = 0; i < n; ++i)
    {
      A1[i] /= A10;
      b[i] /= A10;
    }

  y[0] = -r[0];
  x[0] =  b[0];
  alpha = y[0];
  mult = 1.0 - alpha*alpha;
  beta = 1.0;
  for (size_t k = 1; k < n; ++k)
    {
      size_t k1 = k-1;
      size_t k2 = (k+1)/2 - k%2;

      beta *= mult;
      mu = (b[k] - reverse_inner(k, r, x))
        / beta;
      for (size_t i = 0; i < k; ++i)
          x[i] += mu * y[k1-i];
      x[k] = mu;
      if (k < n1)
        { /* Durbin recursion */
          alpha =
            - (r[k] + reverse_inner(k, r, y))
            / beta;
          mult = 1.0 - alpha*alpha;
          for (size_t i = 0; i < k2; ++i)
            {
              y[i] += alpha * y[k1-i];
              y[k1-i] *= mult;
              y[k1-i] += alpha * y[i];
            }
          if (k%2 == 1)
              y[k2] *= (1.0 + alpha);
          y[k] = alpha;
        }        
    }
  free(y);

  for (size_t i = 0; i < n; ++i)
    {
      A1[i] *= A10;
      b[i] *= A10;
    }

  return 0;
}
