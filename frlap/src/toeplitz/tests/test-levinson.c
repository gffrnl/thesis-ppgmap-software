#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "levinson.h"

#define MAXERR 1.0e-8

int
static test (n, x, y, maxerr)
    const size_t n;
    const double x[restrict const static n];
    const double y[restrict const static n];
    const double maxerr;
{
  double err = 0;

  for (size_t i = 0; i < 3; i++)
    {
      double abserr = fabs(x[i] - y[i]);
      //(void) printf("abserr for i=%ld is %.5f\n", i, abserr);
      err = (abserr > err) ? abserr : err;
    }

  (void) printf("absolute error is %.5e.\n", err);

  if (err > fabs(maxerr))
    {
      (void) puts("Test NOT passed.");
      return EXIT_FAILURE;
    }

  (void) puts("Test passed."); 
  return EXIT_SUCCESS;
}

void
print_arr(n, x)
    const size_t n;
    const double x[const static n];
{
  (void) printf("[");
  for (size_t i = 0; i < n; i++)
      (void) printf("%.5f ", x[i]);
  (void) printf("\b]\n");
  fflush(stdout);
}

int
main (void)
{
  double xe[3] = {355.0, -376.0, 285.0};
  double A1[3] = {10.0, 5.0, 2.0};
  double b[3] = {2240.0, -560.0, 1680.0};
  double x[3];
  
  if (levinson_in_place(3, A1, x, b))
    {
      (void) fprintf(stderr, "Error: null diagonal\n");
      return EXIT_FAILURE;
    }

  (void) printf("A1 = "); print_arr(3, A1);
  (void) printf("b  = "); print_arr(3, b);
  (void) printf("x  = "); print_arr(3, x);

  return test(3, x, xe, MAXERR);
}
