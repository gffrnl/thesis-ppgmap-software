/* frlap/frlap1gd.c */
#include <frlap/frlap1gd.h>
#include <libtoep/toep.c>
#include <assert.h>

static
void
matvecmul_sum (const size_t m,
               const size_t n,
               const double A[const static m][n],
               const double x[const static n],
                     double b[const static m])
{
  for (size_t i = 0; i < m ; i++)
    for (size_t j = 0; j < n; j++)
      b[i] += A[i][j] * x[j];
}

void
tfrlap1gd0(size_t n,
           const double y [static const n],
           size_t ja, size_t jb,
           size_t nc, // must be FRLAP1DQ0_DIFF_COEFF_SIZE(n, ja, jb)
           const double mu[static const nc],
           size_t n0, // must be FRLAP1DQ0_FRLAP_SIZE(n, ja, jb)
           double FLY[static const n0])
{
  /*
   * REMARK:
   *   The products Ba*Ya and Bb*Yb are summed to
   *   FLY before A*Yint for numerical reasons,
   *   because it is expected the values in Ba*Ya
   *   and Bb*Yb are smaller than those in A*Yint
   *   due the decay of the weights.
   */

  size_t na, nb;

  assert(ja > 0 );
  assert(jb < n );
  assert(ja <= jb);
  assert(nc == FRLAP1GD0_DIFF_COEFF_SIZE(n, ja, jb));
  assert(n0 == FRLAP1GD0_FRLAP_SIZE(n, ja, jb));

  na = ja;
  nb = n-1-jb;
  
  {
    /*
     * 4. Assembly of matrices Ba and Bb,
     *    computation of Ba*Ya e Bb*Yb
     *    and sum them to FLY
     */
  
    double (*Ba)[na] = NULL; // to store matrix Ba
    double (*Bb)[nb] = NULL; // to store matrix Bb
    double  *Ya      = NULL; // to store vector Ya
    double  *Yb      = NULL; // to store vector Yb

    // 4.1. Allocation of Ba:
    Ba = (double (*)[na]) malloc(n0 * sizeof(double[na]));
    assert(Ba != NULL); // TODO: change to an error
    // 4.2. Assembly of Ba:
    for (size_t i = 0; i < n0; i++)
      for (ptrdiff_t j = 0; j < na; j++)
        Ba[i][j]=mu[ja+i-j];
    // 4.3. Allocation of Ya:
    Ya = (double *) malloc(na * sizeof(double));
    assert(Ya != NULL); // TODO: change to an error
    // 4.4. Assembly of Ya:
    for (size_t i = 0; i < na; i++)
      Ya[i] = y[i];
    // 4.3. Product Ba*Ya and sum to FLY:  
    matvecmul_sum(n0, na, (const double (*)[na]) Ba, Ya, FLY);
    // 4.4. Deallocation of Ya and Ba:
    free(Ya); // *** array Ya no more needed ***
    free(Ba); // *** array Ba no more needed ***
    
    // 4.5. Allocation of Bb:
    Bb = (double (*)[nb]) malloc(n0 * sizeof(double[nb]));
    assert(Bb != NULL); // TODO: change to an error
    // 4.6. Assembly of Bb:
    for (ptrdiff_t i = 0; i < n0; i++)
      for (size_t j = 0; j < nb; j++)
        Bb[i][j]=mu[n0-i+j];
    // 4.3. Allocation of Yb:
    Yb = (double *) malloc(nb * sizeof(double));
    assert(Yb != NULL); // TODO: change to an error
    // 4.4. Assembly of Yb:
    for (size_t i = 0; i < nb; i++)
      Yb[i] = y[i+jb+1];
    // 4.7. Product Bb*Yb and sum to FLY:
    matvecmul_sum(n0, nb, (const double (*)[nb]) Bb, Yb, FLY);
    // 4.8. Deallocation of Yb and Bb:
    free(Yb); // *** array Yb no more needed ***
    free(Bb); // *** array Bb no more needed ***
  }

  {
    /*
     * 5. Assembly of matrices Yint, computation
     *    of A*Yint and sum it to FLY
     */

    double *Yint  = NULL; // to store vector Yint
    double *AYint = NULL; // To store vector A*Yint

    // 5.1. Allocation of Yint:
    Yint = (double *) malloc(n0 * sizeof(double));
    assert(Yint != NULL); // TODO: change to an error
    // 5.1. Assembly of Yint:
    for (size_t i = 0; i < n0; i++)
      Yint[i] = y[i+ja];
    // 5.3. Allocation of AYint:
    AYint = (double *) malloc(n0 * sizeof(double));
    assert(AYint != NULL); // TODO: change to an error
    // 5.4. Fast symmetric toeplitz-vector A*Yint:
    int ret = fast_symm_toeplitz_prod(n0, mu, Yint, AYint);
    assert(ret == 0); // TODO: change to an error
    // 5.5. Deallocation of Yint:
    free(Yint); // *** array Yint no more needed ***
    // 5.6. Summation of AYint to FLY:
    for (size_t i = 0; i < n0; i++)
      FLY[i] += AYint[i];
    // 5.7. Deallocation of AYint:
    free(AYint); // *** array AYint no more needed ***
  }
}
