#include <stddef.h>

int
toeplitz_prod (size_t n,
    const double A[const restrict static n],
    const double x[const restrict static n],
          double b[const restrict static n]);
{
  //*
  for (size_t i = 0; i < n; i++)
    {
      h[i] = 0.0;
      for (size_t j = n-(i+1); j > 0; j--)
          h[i] += f[j] * g[j+i];
      for (size_t j = i-1; j < i; j--)
          h[i] += f[i-j] * g[j];
      h[i] += f[0] * g[i];
    }
  //*/
  /*
  size_t i;
  for (i = 0; i < n-(i+1); i++)
    {
      h[i] = 0.0;
      for (size_t j = n-(i+1); j > i; j--)
          h[i] += f[j] * g[i+j];
      for (size_t j = i; j > 0; j--)
          h[i] += f[j] * (g[i-j] + g[i+j]);
      h[i] += f[0] * g[i];
      
    }
  for (; i < n; i++)
    {
      h[i] = 0.0;
      for (size_t j = i; j > n-(i+1); j--)
      {
          h[i] += f[j] * g[i-j];
      }
      for (size_t j = n-(i+1); j > 0; j--)
          h[i] += f[j] * (g[i-j] + g[i+j]);
      h[i] += f[0] * g[i];
    }
  //*/
  
  return 0;
}

