/* frlap/frlap1dq/wgtor_centered_periodic_3_point.c */
//#include "wgtor.h" /* maybe change to <frlap/frlap1qd/wgtor/wgtor.h> */
#include <math.h>
#include <stddef.h>

/*
 * TODO:
 *   - maybe put the beta function in a library of
 *     special math functions
 */

/**
 *  Math special function beta(a, b);
 *    - here a, b are reals, but in fact can be complex
 *    - definition:
 *        beta(a, b) =
 *          (gamma(a) * gamma(b)) / gamma(a + b) 
 */
static double beta (double a, double b)
{
 /*
  * REMARK: using log-gamma function and then exponentiate
  *         is better to mitigate round-off/approximation
  *         errors
  */
  return exp(lgamma(a) + lgamma(b) - lgamma(a + b));
}

void
frlap1dq_wgtor_centered_periodic_3_point (double order,
                                          double gstep,
                                          size_t n,
                                          double w[const static n])
{
  double c1 = 0.0, c2 = 0.0;

  c1 = 1.0 / (pow(gstep, - order) * M_PI);
  c2 = c1 * sin(0.5 * order * M_PI);
  
  w[0] = - c1 * exp2(order) * beta(0.5 + 0.5 * order, 0.5);

  for (size_t k = 1; k < n; k++)
    w[k] = c2 * beta((double) k - 0.5 * order, 1.0 + order);
}
