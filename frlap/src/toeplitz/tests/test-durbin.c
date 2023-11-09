#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "durbin.h"

#define MAXERR 1.0e-8
#ifndef N
  #define N 3
#endif

int main(int argc, char *argv[])
{
  double r[] = {0.5, 0.2, 0.1, 0.9, 0.3};
  double x[N]
#if N == 3
    = { -75.0/140,
         12.0/140,
        - 5.0/140};
#elif N == 4
    = { - 517.0/1044,
        -  12.0/1044,
          597.0/1044,
        -1184.0/1044};
#elif N == 5
    = {  6258780.0/2326032,
        -3763359.0/2326032,
         1405224.0/2326032,
          597951.0/2326032,
        -6534396.0/2326032 };
#else
  ;
  #error minimum value for N is 3 and maximum is 5
#endif
  double y[N];
  double err;
  
  durbin(N, r, y);
  err=0;
  for (size_t i = 0; i < N; i++)
    {
      double abserr = fabs(x[i] - y[i]);
      //(void) printf("abserr for i=%ld is %.5f\n", i, abserr);
      err = (abserr > err) ? abserr : err;
    }
  (void) printf("absolute error is %.5e.\n", err);
  if (!(err > MAXERR))
    {
      (void) puts("Test passed.");
      return EXIT_FAILURE;
    }
  else
      (void) puts("Test NOT passed."); 
  
  return EXIT_SUCCESS;
}
