#pragma once

#include <cstddef>
#include <cstring>
#include <cmath>
#include <vector>
#include <fftw3.h>

std::vector<double>
fast_symm_prod(std::vector<double> const & c,
               std::vector<double> const & x) {

  // TODO: VALIDATE ARGUMENTS

  std::size_t n = c.size();
  std::size_t k; // padding
  std::size_t m; // the dimension of the augmented matrix

  // Augmentation:
  // A1 and x must be embedded into n+2*(k+1) arrays
  std::vector<double> mu; // first row of augmented matrix
  std::vector<double> y;  // augmented vector
  std::vector<double> aux; // auxiliary array

  fftw_plan plan;

  // The sizes needed to allocate memory for the augmented arrays
  k = (n % 2 == 0) ? (n - 4) / 2 : (n - 3) / 2;
  m = n + 2 * (k + 1);

  // Augmented arrays allocation
  mu.resize(m, 0);
  y.resize(m, 0);
  aux.resize(m, 0);

  // Construct the augmented arrays
  std::memcpy((void *) mu.data(), (const void *) c.data(), n * sizeof(double));
  for (size_t i = 2; i < n; ++i)
      mu[i-2] -= c[i];
  std::memcpy((void *) (y.data()+k+1), (const void *) x.data(), n * sizeof(double));

  // Compute the DST1 of mu and store in aux
  plan =
    fftw_plan_r2r_1d(m, mu.data(), aux.data(), FFTW_RODFT00, FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);
  fftw_cleanup();

  // Compute the DST1 of y and store in mu
  plan = fftw_plan_r2r_1d(m, y.data(), mu.data(), FFTW_RODFT00, FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);
  fftw_cleanup();

  y.resize(0);
  y.shrink_to_fit(); // y no more needed

  // Multiply
  {
    double scaling = 1.0 / ( 4.0 * (m + 1) );
    for (std::size_t i = 0; i < m; ++i)
      aux[i] *= (scaling * mu[i] / std::sin( (i+1) * M_PI / (m+1) ));
  }
  // NOTE: need to multiply lamb by 1/2 rather
  //         sqrt(n+1)/2 because in fftw3 scaling
  //         of DST1 is 2
  
  // Now compute the DST1 of xemb
  plan =
    fftw_plan_r2r_1d(m, aux.data(), mu.data(), FFTW_RODFT00, FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);
  fftw_cleanup();

  aux.resize(0);
  aux.shrink_to_fit(); // aux no more needed

  std::vector<double> b(n);
  
  /* Copy only necessary values to b */
  std::memcpy((void *) b.data(), (const void *) (mu.data()+k+1), n * sizeof(double));

  return b;
}
