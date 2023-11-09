/* test of intde1.c */

#include <math.h>
#include <stdio.h>

#include "intde1.c"

#ifndef N
  #define N 5
#endif

#define K 10

int k = 0;

int main()
{
  extern int k;

  double f(double);
  void intde(double (*f)(double), double a, double b, double eps, double *i, double *err);
  double ei, i, err, cerr;

  for (int l = 0; l < K; l++)
    {
      k = l;
      if (k == 0)
          ei = M_PI*M_PI*M_PI/3.0;
      else if (k%2 == 0)
          ei =   (2*M_PI) / (double) (k*k);
      else
          ei = - (2*M_PI) / (double) (k*k);
      intde(f, 0.0, M_PI, 1.0e-15, &i, &err);
      cerr = ei - i;
      printf("%d \t %+.5e \t %+.5e \t %+.5e \t %+.5e\n",
             k, ei, i, err, cerr);
    }

  return 0;
}


double f(double x)
{
  extern int k;
  
  return x*x * cos(k*x);
}

