/* frlap/frlap1dq/wgtor/wgtor/huang_oberman_linear.c */
//#include "wgtor.h"
#include <stddef.h>
#include <math.h>
#include <float.h>

#include <stdbool.h> // TODO: Maybe remove

#include <stdio.h> // TODO: Maybe remove

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

static double frlap_ncons(size_t n, double order)
{
  return order * exp2(order - 1)
    / pow(M_PI, 0.5 * n)
      * exp(  lgamma(0.5 * (order + n))
            - lgamma( 1.0 - 0.5*order )  );
}

static double frlap_ncons1(double order)
{
  return frlap_ncons(1, order);
}

void
frlap1dq_wgtor_huang_oberman_linear (double order,
                                     double gstep,
                                     size_t n,
                                     double w[const static n])
{
  extern double d1G_alpha_ne_1 (double, size_t);
  extern double d2G_alpha_ne_1 (double, size_t);

  extern double d1G_alpha_eq_1 (size_t);
  extern double d2G_alpha_eq_1 (size_t);

  double c1 = 0.0, c2 = 0.0;

  puts("\n ** in function frlap1qd_wgtor_huang_oberman_linear **");
  fflush(stdout);

  if (almost_equal(order, 1.0))
    c1 = 1.0 / gstep;
  else
    c1 = pow(gstep, - order);

  c2 = frlap_ncons1(order) * c1;
  
 /*
  * REMARK: using log-gamma function and then
  *         exponentiate is better to mitigate
  *         round-off/approximation errors
  */
  w[0] = - c1 * exp2(order) / sqrt(M_PI)
    * exp(  lgamma(0.5 * (1.0 + order))
          - lgamma(2.0 - 0.5 * order  ) );

  if (almost_equal(order, 1.0))
    {
      w[1] = c2 * ( 1.0 - d2G_alpha_eq_1(1)
                        + d1G_alpha_eq_1(2)
                        - d1G_alpha_eq_1(1) );

      for (size_t k = 2; k < n; k++)
        w[k] =
          c2 * (         d1G_alpha_eq_1(k+1)
                 +       d1G_alpha_eq_1(k-1)
                 - 2.0 * d1G_alpha_eq_1(k)   );
    }
  else
    {
      w[1] = c2 * (   1.0 / (2.0 - order)
                    - d2G_alpha_ne_1(order, 1)
                    + d1G_alpha_ne_1(order, 2)
                    - d1G_alpha_ne_1(order, 1) );

      for (size_t k = 2; k < n; k++)
        w[k] =
          c2 * (         d1G_alpha_ne_1(order, k+1)
                 +       d1G_alpha_ne_1(order, k-1)
                 - 2.0 * d1G_alpha_ne_1(order, k)   );
    }
}
