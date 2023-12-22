#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>

#include <frlap/frlap1gd.h>

#include <array.h>
#include <linspace.h>
#include <fltcmp.h>

#include <assert.h>

/*
#include "array.h"
#include "linspace.h"
#include "fltcmp.h"
#include "frlap/frlap1dq.h"
*/

/*
#include <assert.h>
#include <time.h>
*/

/*
double
maxabserr(const size_t n,
          const double x[const restrict static n],
          const double y[const restrict static n],
          const size_t imin,
          const size_t imax,
                size_t *k)
{
  assert(imin >= 0);
  assert(imax <  n);

  double abserr = 0.0;

  for (size_t i = imin; i <= imax; i++)
    {
      double temp = fabs(x[i] - y[i]);
      if (temp > abserr)
        {
          abserr = temp;
          *k = i;
        }
    }
  return abserr;
}
*/

/*
double
maxrelerr(const size_t n,
          const double x[const restrict static n],
          const double y[const restrict static n],
          const size_t imin,
          const size_t imax,
                size_t *k)
{
  assert(imin >= 0);
  assert(imax <  n);

  double relerr = 0.0;

  for (size_t i = imin; i <= imax; i++)
    {
      double temp = fabs( (x[i] - y[i]) / fmin(x[i], y[i]) );
      if (temp > relerr)
        {
          relerr = temp;
          *k = i;
        }
    }
  return relerr;
}
*/

static double alpha = 0.4;

int
main (void)
{
  double f (double); // the function we want to compute the frac. Laplacian
  double frLap_f (double); // exact frac. Laplacian of function f
  double a, b; // the endpoints of the computational grid
  double h; // the grid spacing
  size_t n; // number of grid points
  
  double *x, *y, *frLapy;

  size_t n0, ja, jb;
  double xja, xjb;

  //struct frlap1gd prm;
  
  // ----VERIFICAR----
  /*
  double *x, *y, *frLapfx, *frLapy, *frLapy_naive;
  struct frlap1dq prm;
  clock_t t;
  */
  // ----VERIFICAR----

  b = 6.0;
  a = -b;
  n = 121;

  xjb = 2.0;
  xja = -xjb;

  // construct the grid
  //x = (double *) malloc(n * sizeof(double));
  x = (double *) malloc(n * sizeof(double));
  if (x == NULL)
    {
      perror("malloc() failed for x");
      exit(1);
    }
  h = linspace(a, b, n, x);
  assert(dblcmp_dig(h, 0.1, 6)); // asserts if h is equal 0.1 up to 6 digits os precision
  
  /*
  (void) printf("grid spacing h is %f\n", h);
  (void) puts("x:");
  print_arr(n, x);
  (void)puts(" ");
  (void) fflush(stdout);
  */

  {
    ptrdiff_t sja, sjb;

    sja = seekval_dig(xja, n, x, 6);
    if (sja == -1)
      {
        (void) fprintf(stderr, "cannot find xja=%f in x\n", xja);
        exit(3);
      }
    else
      ja = (size_t) sja;

    sjb = seekval_dig(xjb, n, x, 6);
    if (sja == -1)
      {
        (void) fprintf(stderr, "cannot find xjb=%f in x\n", xjb);
        exit(3);
      }
    else
      jb = (size_t) sjb;

    n0 = jb - ja + 1;
    (void) printf("ja = %lld, x[ja] = %f\n", ja, xja);
    (void) printf("jb = %lld, x[jb] = %f\n", jb, xjb);
    (void) printf("n0 = %lld\n", n0);
  }

  y = (double *) malloc(n * sizeof(double));
  if (y == NULL)
    {
      perror("malloc() failed for y");
      free(x);
      exit(2);
    }
  elemental(n, x, f, y);

  
  /*
  (void) puts("y:");
  print_arr(n, y);
  (void)puts(" ");
  (void) fflush(stdout);
  */

  /*
  prm.order = alpha;
  prm.grid_step = h;
  prm.coeff_type = FRLAP1GD_COEFF_TYPE_HUANG_OBERMAN_QUADRATIC;

  tfrlap1gd(prm, n, y, frLapy);
  */
  /*
  for (unsigned expon = 6; expon < 25; expon++)
    {
      double a, b;
      size_t kabs, krel;

      n = (1 << expon) - 1; // 2^3 = 8

      (void) printf("n = %lld\n", n);
      fflush(stdout);

      a = - ((double) (n-1)) / 20;
      b =   ((double) (n-1)) / 20;

      (void) printf("[a, b] = [%f, %f]\n", a, b);
      fflush(stdout);
      assert(isless(a, -2.0) && isgreater(b, 2.0));

      x = (double *) malloc(n * sizeof(double));
      y = (double *) malloc(n * sizeof(double));
      frLapfx =
        (double *) malloc(n * sizeof(double));
      frLapy =
        (double *) malloc(n * sizeof(double));
      frLapy_naive =
        (double *) malloc(n * sizeof(double));

      assert(x            != NULL);
      assert(y            != NULL);
      assert(frLapfx      != NULL);
      assert(frLapy       != NULL);
      assert(frLapy_naive != NULL);
      
      h = linspace(a, b, n, x);
      assert(dblcmp_dig(h, 0.1, 6) == true);

      (void) printf("x[n/2-21] = %f\n", x[n/2-21]);
      (void) printf("x[n/2+20] = %f\n", x[n/2+20]);
      fflush(stdout);
      
      elemental(n, x, f, y);
      elemental(n, x, frLap_f, frLapfx);

      prm.order = alpha;
      prm.gstep = h;
      prm.wtype =
        FRLAP1DQ_WTYPE_HUANG_OBERMAN_QUADRATIC;

      t = clock();
      //tfrlap1dq__naive__(prm, n, y, frLapy);
      tfrlap1dq(prm, n, y, frLapy);
      t = clock() - t;

      
      (void) printf("\ttime to compute the truncated fractional Laplacian: ");
      (void) printf("%f s\n", ((double) t)/CLOCKS_PER_SEC);
      
      (void) printf("\tmaximum absolute error: %.5e ",
                    maxabserr(n, frLapfx, frLapy, n/2-21, n/2+20, &kabs));
      (void) printf(" (at x = %f)\n", x[kabs]);
      (void) printf("\tmaximum relative error: %.5e ",
                    maxrelerr(n, frLapfx, frLapy, n/2-21, n/2+20, &krel));
      (void) printf(" (at x = %f)\n", x[krel]);

      (void) puts("------");
      fflush(stdout);
  */
  //{ /* Save to file */
        /*
        FILE *fp;
        fp = fopen("bench.dat", "w");
        assert(fp != NULL);
        fprintf(fp, "#%16s\t%16s\t%16s\t%16s\t%16s\t%16s\t\n",
                "x", "f(x)", "frLap-approx", "frLap-exact", "abserr", "relerr");
        for (size_t i = 0; i < n; i++)
          fprintf(fp, "%16.8e\t%16.8e\t%16.8e\t%16.8e\t%16.8e\t%16.8e\t\n",
                  x[i], y[i], frLapy[i], frLapfx[i],
                  fabs(frLapy[i]-frLapfx[i]),
                  fabs((frLapy[i]-frLapfx[i])/frLapfx[i]));
        fclose(fp);
        */
  // }

      /*
      free(frLapy_naive);
      free(frLapy);
      free(frLapfx);
      free(y);
      free(x);
      */
      //}

  free(y);
  free(x);
  return 0;
}

double
f (double x)
{
  extern double alpha;
  return pow(1.0 + x*x, -(1.0-alpha)/2.0);
}

double
frLap_f (double x)
{
  extern double alpha;
  return
    exp2(alpha) *
    exp(  lgamma((1.0+alpha)/2.0)
        - lgamma((1.0-alpha)/2.0) ) *
    pow(1.0 + x*x, -(1.0+alpha)/2.0);
}
