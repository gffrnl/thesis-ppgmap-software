#include <stddef.h>

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
