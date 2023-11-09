#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>

#include "array.h"
#include "linspace.h"
#include "fltcmp.h"
#include <frlap/frlap1dq.h>

#include <assert.h>


static double alpha = 0.4;

int
main (void)
{
  double f (double);
  double frLap_f (double);

  double a, b;
  size_t n;

  double a0, b0;
  size_t ja0, jb0;

  double h, *x, *y, *frLapfx, *frLapy;
  struct frlap1dq prm;

  a = -6.0;
  b =  6.0;
  n = 121;
  
  a0 = -2.0;
  b0 =  2.0;


  
  x = (double *) malloc(n * sizeof(double)); // TODO: check if x is allocated!
  y = (double *) malloc(n * sizeof(double)); // TODO: check if y is allocated!
  frLapfx =
    (double *) malloc(n * sizeof(double)); // TODO: check if z is allocated!  
  frLapy =
    (double *) malloc(n * sizeof(double)); // TODO: check if it is allocated!

  h = linspace(a, b, n, x);
  assert(dblcmp_dig(h, 0.1, 6) == true);

  elemental(n, x, f, y);
  elemental(n, x, frLap_f, frLapfx);

  prm.order = alpha;
  prm.gstep = h;
  prm.wtype =
    FRLAP1DQ_WTYPE_HUANG_OBERMAN_QUADRATIC;

  assert(a0 >= a);
  assert(b0 <= b);

  ja0 = (a0 - a) / h;
  assert(dblcmp_dig(a + ja0*h, a0, 6) == true);
  jb0 = (b0 - a) / h;
  assert(dblcmp_dig(a + jb0*h, b0, 6) == true);

  //tfrlap1dq(prm, n, y, frLapy);
  tfrlap1dq_0(prm, n, y, ja0, jb0, frLapy+ja0);
  
  { /* Save to file */
    FILE *fp;
    fp = fopen("bench_0.dat", "w");
    assert(fp != NULL);
    fprintf(fp, "#%16s\t%16s\t%16s\t%16s\t%16s\t%16s\t\n",
            "x", "f(x)", "frLap-approx", "frLap-exact", "abserr", "relerr");
    for (size_t i = ja0; i <= jb0; i++)
      fprintf(fp, "%16.8e\t%16.8e\t%16.8e\t%16.8e\t%16.8e\t%16.8e\t\n",
              x[i], y[i], frLapy[i], frLapfx[i],
              fabs(frLapy[i]-frLapfx[i]),
              fabs((frLapy[i]-frLapfx[i])/frLapfx[i]));
    fclose(fp);
  }

  free(frLapy);
  free(frLapfx);
  free(y);
  free(x);

  return 0;
}

double
f (double x)
{
  extern double alpha;
  return pow(1.0 + x*x, alpha/2.0-0.5);
}

double
frLap_f (double x)
{
  extern double alpha;
  return
      pow(2.0, alpha) *
      exp(  lgamma(0.5+alpha/2.0)
          - lgamma(0.5-alpha/2.0) ) *
      pow(1.0 + x*x, -0.5-alpha/2.0);
}
