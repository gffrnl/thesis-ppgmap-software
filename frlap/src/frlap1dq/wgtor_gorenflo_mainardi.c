/* frlap/frlap1dq/wgtor_gorenflo_mainardi.c */
//#include "wgtor.h"
#include <math.h>
#include <stddef.h>
#include <stdbool.h>
#include <float.h>

static bool almost_equal(double a, double b)
{
  return fabs(a - b) < DBL_EPSILON;
}


void
frlap1dq_wgtor_gorenflo_mainardi (double order,
                                  double gstep,
                                  size_t n,
                                  double w[const static n])
{
  /*
   * TODO:
   *   - seek for round-off errors, catastrophic
   *     cancelation or other numeric problems
   *   - maybe rewrite for optimization
   */

  double c1 = 0.0, c2 = 0.0;

  if (almost_equal(order, 1.0))
    c1 = 1.0 / (gstep * M_PI);
  else
    {
      c1 = - pow(gstep, - order)
        / ( 2.0 * cos(0.5 * order * M_PI) );

      c2 = c1 * order * (order - 1)
        / tgamma(2.0 - order);
    }

  if (almost_equal(order, 1.0))
    {
      w[0] = - 2.0 * c1;

      for (size_t k = 1; k < n; k++)
          w[k] = c1 / (k * (k+1));
    }
  else if (order < 1.0)
    {
      w[0] = 2.0 * c1;

      /*
       * REMARK: using log-gamma function and then
       *         exponentiate is better to mitigate
       *         round-off/approximation errors
       */
      for (size_t k = 1; k < n; k++)
        w[k] = c2 * exp( lgamma(k-order) - lgamma(k+1) );
    }
  else if (order > 1.0)
    {
      w[0] = - 2.0 * order * c1;

      w[1] = c1 * (1.0 + 0.5 * order * (order - 1.0));

      /*
       * REMARK: using log-gamma function and then
       *         exponentiate is better to mitigate
       *         round-off/approximation errors
       */
      for (size_t k = 2; k < n; k++)
        w[k] = c2 * exp( lgamma(k+1-order) - lgamma(k+2) );
    }
}
