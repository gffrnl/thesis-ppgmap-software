/* frlap/frlap1dq/wgtor_huang_oberman.h */
#include <stddef.h>
#include <math.h>

double
d0G_alpha_ne_1 (double alpha, size_t k)
{
  return - pow((double) k, 2.0 - alpha)
    / ( (2.0 - alpha) * (1.0 - alpha) * alpha );
}

double
d1G_alpha_ne_1 (double alpha, size_t k)
{
  return - pow((double) k, 1.0 - alpha)
    / ( (1.0 - alpha) * alpha );
}

double
d2G_alpha_ne_1 (double alpha, size_t k)
{
  return - pow((double) k, -alpha) / alpha;
}

double
d0G_alpha_eq_1 (size_t k)
{
  return (double) k - k * log((double) k);
}

double
d1G_alpha_eq_1 (size_t k)
{
  return - log((double) k);
}

double
d2G_alpha_eq_1 (size_t k)
{
  return - 1.0 / (double) k;
}
