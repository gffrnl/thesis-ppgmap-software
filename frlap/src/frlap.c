#include <frlap.h>
#include <math.h>

double
frlap_ncons (size_t n, double order)
{
  return order * exp2(order - 1)
    / pow(M_PI, 0.5 * n)
      * exp(  lgamma(0.5 * (order + n))
            - lgamma( 1.0 - 0.5*order )  );
}

double frlap_ncons1 (double order)
{
  return frlap_ncons(1, order);
}
