/* frlap/frlap1dq/wgtor_periodic_centered_5_point.c */
//#include "wgtor.h" /* maybe change to <frlap/frlap1qd/wgtor/wgtor.h> */
//#include <math.h>
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

//static double beta (double a, double b)
//{
// /*
//  * REMARK: using log-gamma function and then exponentiate
//  *         is better to mitigate round-off/approximation
//  *         errors
//  */
//  return exp(lgamma(a) + lgamma(b) - lgamma(a + b));
//}

#include <stdio.h>
#include <stdlib.h>
void
frlap1dq_wgtor_centered_periodic_5_point (double order,
                                          double gstep,
                                          size_t n,
                                          double mu[const static n])
{
  (void) fprintf(stderr,
    "frlap1qd_wgtor_periodic_centered_5_point() not implemented\n");
  abort();
}
