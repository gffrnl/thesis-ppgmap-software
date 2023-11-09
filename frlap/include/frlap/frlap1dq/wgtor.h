/* frlap/frlap1dq/wgtor.h */
#pragma once
#include <stddef.h>

void
frlap1dq_wgtor_spectral (double order,
                         double gstep,
                         size_t n,
                         double w[const static n]);
void
frlap1dq_wgtor_huang_oberman_linear (double order,
                                     double gstep,
                                     size_t n,
                                     double w[const static n]);
void
frlap1dq_wgtor_huang_oberman_quadratic (double order,
                                        double gstep,
                                        size_t n,
                                        double w[const static n]);

void
frlap1dq_wgtor_gorenflo_mainardi (double order,
                                  double gstep,
                                  size_t n,
                                  double w[const static n]);

void
frlap1dq_wgtor_centered_periodic_3_point (double order,
                                          double gstep,
                                          size_t n,
                                          double w[const static n]);

void
frlap1dq_wgtor_centered_periodic_5_point (double order,
                                          double gstep,
                                          size_t n,
                                          double w[const static n]);

/*
void
muper3_ooura (
    const double alpha,
    const double dx,
    const size_t n,
          double mu[const static n],double err[const static n]);
*/          
