#ifndef PSKE_H
#define PSKE_H
#ifndef __cplusplus
#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <stdbool.h> 
#endif
#ifdef __cplusplus
extern "C" {
#endif
  
/* declare functions */
// functions in skeletoid_expm_dense.c
SEXP R_skeletoid_expm_dense(SEXP Q_, SEXP delta_, SEXP n_, SEXP n_sq_);
void skeletoid_expm_dense(double *S, double *Q, double delta, size_t n, size_t n_sq);
void build_X_delta(double *X, double *Q, double delta, size_t n);
void matsq_hp(int n, int k, double *X);
// void matsq(int n, int k, double *z);

// functions in skeletoid_vtexpm_dense.c
SEXP R_skeletoid_vtexpm_dense(SEXP Q_, SEXP delta_, SEXP v_,
                              SEXP nr_, SEXP nc_, SEXP k1_, SEXP k2_);
void skeletoid_vtexpm_dense(double *A, double *Q, double delta, 
                            double *v, size_t nr, size_t nc, int k1, int k2);
void mat_mat_prod_loop(double *z, double *x, double *v, int nr, int nc, int k,
                       bool do_init);

#ifdef __cplusplus
}
#endif
#endif
