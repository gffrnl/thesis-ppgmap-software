/* toeplitz.c */
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <fftw3.h>

void
symm_toeplitz_prod (const size_t n,
    const double A1[const restrict static n],
    const double x [const restrict static n],
          double b [const restrict static n])
{
  for (size_t i = 0; i < n; ++i)
    {
      b[i] = 0.0;
      for (size_t j = n-(i+1); j > 0; j--)
          b[i] += A1[j] * x[j+i];
      for (size_t j = i-1; j < i; j--)
          b[i] += A1[i-j] * x[j];
      b[i] += A1[0] * x[i];
    }
}

int
fast_symm_toeplitz_prod (const size_t n,
    const double A1[const restrict static n],
    const double x [const restrict static n],
          double b [const restrict static n])
{
  size_t k;   // padding
  size_t m;   /* the dimension of the augmented
                 matrix */
  
  /* Augmented arrays: A1 and x must be embedded
     into n+2*(k+1) arrays */
  double *mu;  // first row of augmented matrix
  double *y;   // augmented vector
  double *aux; // auxiliary array

  fftw_plan plan;

  /* The sizes needed to allocate memory for the
     augmented arrays */
  k = (n % 2 == 0) ? (n - 4) / 2 : (n - 3) / 2;
  m = n + 2 * (k + 1);

  /* Augmented arrays allocation */
  /*
  if ((mu  = (double *) malloc(m * sizeof(double))) == NULL)
      return 1;
  if ((y   = (double *) malloc(m * sizeof(double))) == NULL)
      return 2;

  if ((aux = (double *) malloc(m * sizeof(double))) == NULL)
      return 3;
  */
  
  if ((mu  = (double *) calloc(m, sizeof(double))) == NULL)
      return 1;
  if ((y   = (double *) calloc(m,  sizeof(double))) == NULL)
      return 2;

  if ((aux = (double *) calloc(m, sizeof(double))) == NULL)
      return 3;


  /* Construct of the the augmented arrays */
  /*
  memset((void *) (mu+n)   , 0, (m-n) * sizeof(double));
  memset((void *)  y       , 0, (k+1) * sizeof(double));
  memset((void *) (x+k+n+1), 0, (k+1) * sizeof(double));
  // TODO: INSTEAD THE ABOVE, TRY CALLOC
  // WARNING: THESE memsets are wrong !!!
  */
  
  memcpy((void *) mu, (const void *) A1, n * sizeof(double));
  for (size_t i = 2; i < n; ++i)
      mu[i-2] -= A1[i];
  memcpy((void *) (y+k+1), (const void *) x, n * sizeof(double));

  /* Compute the DST1 of mu and store in aux */
  plan = fftw_plan_r2r_1d(m, mu, aux,
    FFTW_RODFT00, FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);
  fftw_cleanup();

  /* Compute the DST1 of y and store in mu */
  plan = fftw_plan_r2r_1d(m, y, mu,
    FFTW_RODFT00, FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);
  fftw_cleanup();

  free(y); // y no more needed

  /* Multiply */
  {
    double scaling = 1.0 / ( 4.0 * (m + 1) );
    for (size_t i = 0; i < m; ++i)
        aux[i] *= (scaling * mu[i] / sin( (i+1) * M_PI / (m+1) ));
  }
  /* NOTE: need to multiply lamb by 1/2 rather
           sqrt(n+1)/2 because in fftw3 scaling
           of DST1 is 2 */

  /* Now compute the DST1 of xemb */
  plan = fftw_plan_r2r_1d(m, aux, mu,
    FFTW_RODFT00, FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);
  fftw_cleanup();

  free(aux); // aux no more needed

  /* Copy only necessary values to b */
  memcpy((void *) b, (const void *) (mu+k+1), n * sizeof(double));

  free(mu);

  return 0;
}

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

void
durbin (n, r, y)
    const size_t n;
    const double r[restrict const static n];
          double y[restrict const static n];
{
  double alpha, beta, mult;

  y[0] = -r[0];
  alpha = y[0];
  mult = 1.0 - alpha*alpha;
  beta = 1.0;
  for (size_t k = 1; k < n; ++k)
    {
      size_t k1 = k-1;
      size_t k2 = (k+1)/2 - k%2;

      beta *= mult;
      alpha = - (r[k] + reverse_inner(k, r, y))
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
