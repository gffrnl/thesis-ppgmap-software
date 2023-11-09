/* frlap/frlap1dq/wgtor_spectral.c */
//#include "wgtor.h"
#include <stddef.h>
//#include "ooura/intde/intde1.c"
#include <stdbool.h>
#include <float.h>
#include <math.h>

/*
 * TODO:
 *   - maybe put the function almost equal in a library
 *   - see  https://en.cppreference.com/w/cpp/types/numeric_limits/epsilon
 *     for a better implementation of the equality between floating point
 *     types
 */

static bool almost_equal(double a, double b)
{
  return fabs(a - b) < DBL_EPSILON;
}

static int freq;
static double expo;

static double f(double x)
{
  extern int freq;
  extern double expo;

  return pow(x, expo) * cos(freq * x);
}

/*
void
frlap1dq_wgtor_spectral (double order,
                         double gstep,
                         size_t n,
                         double w[const static n])
{
  if (almost_equal(order, 1.0))
    {
      w[0] = -1.0;
      for (size_t k = 1; k < n; k+=2)
        w[k] = 2.0 / (M_PI * k*k * gstep);
      for (size_t k = 2; k < n; k+=2)
        w[k] = 0.0;
    }
  else if (almost_equal(order, 2.0))
    {
      w[0] = - (M_PI / gstep) * (M_PI / gstep) / 3.0;
      for (size_t k = 1; k < n; k+=2)
        w[k] =   2.0 / ( (k*h) * (k*h) );
      for (size_t k = 2; k < n; k+=2)
        w[k] = - 2.0 / ( (k*h) * (k*h) );
    }
  else
    {
      extern int freq;
      extern double expo;
      double err;

      expo = order;
      
      w[0] = - pow(M_PI / gstep, order) / (order + 1.0);

      gstep = pow(gstep, order);

      for (size_t k = 1; k < n; k++)
        {
          freq = k;
          intde(f, 0.0, M_PI, 1.1e-14, &mu[k], &err);
          w[k] *= -1.0 / (M_PI * gstep)
        }
    }
}
*/
