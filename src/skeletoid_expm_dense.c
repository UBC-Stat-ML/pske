#define USE_FC_LEN_T
#include "pske.h"
#ifndef FCONE
# define FCONE
#endif

/*
 *  skeletoid_expm high-precision version, deals with rounding issues when delta ~ 0
 *  Basically same cost as the naive version so we just use this for all 
 *  S = I+(S-I) := I+X0. Then
 *  S^2 = (I + X0)^2 = I + 2X0 + X0^2 = I+X1, with X1:=2X0 + X0^2
 *  By induction: S^(2^k) = I+X_k , with X_{k+1} = 2X_k + X_k^2 and X0 = S-I
 */

/* .Call() this from R */
SEXP R_skeletoid_expm_dense(SEXP Q_, SEXP delta_, SEXP n_, SEXP n_sq_) {
  R_len_t n = asInteger(n_), n_sq = asInteger(n_sq_);
  double delta = asReal(delta_);
  int nprotect = 0;
  SEXP S_sexp; // declare answer
  PROTECT(S_sexp = allocMatrix(REALSXP, n, n)); nprotect++; // allocate+PROTECT
  skeletoid_expm_dense(REAL(S_sexp), REAL(Q_), delta, n, n_sq); // populate S_delta and do squaring
  UNPROTECT(nprotect);
  return(S_sexp); /* unprotect ans and return it */
}

void skeletoid_expm_dense(double *X, double *Q, double delta, size_t n, size_t n_sq) {
  build_X_delta(X,Q,delta,n); // build X_delta = S_delta - I
  if(n_sq > 0) // modified squaring loop for X
    matsq_hp(n, n_sq, X);
  for(int i=0; i<n; i++) // Recover S^(2^n_sq) = Id + X^(2^n_sq)
    X[i+n*i]+=1.0;
}

void build_X_delta(double *X, double *Q, double delta, size_t n){
  size_t i, j;
  double dqj, qii_m_qjj;
  double* dQ = (double*)malloc(sizeof(double) * n); // dQ = diag(Q)
  double* p_0 = (double*)malloc(sizeof(double) * n); // p_0 = exp(delta*diag(Q))
  double* S_sym = (double*)malloc(sizeof(double) * n * n); // symmetric part of S_delta
  
  /*  fetch diag(Q), compute diagonal of S_delta=p_0 */
  for( i=0; i<n; i++){
    dQ[i] = Q[i + n*i];
    p_0[i] = exp(delta*dQ[i]);
  }
  
  /* populate upper triangle of S_sym */
  for( j=0; j<n; j++){
    dqj = dQ[j]; // fetch qjj
    for( i=0; i<j; i++){ // only do (strict) upper triangle
      qii_m_qjj = dQ[i] - dqj;
      if(fabs(qii_m_qjj) < 100.0*DBL_EPSILON ) {
        S_sym[i + n*j] = p_0[i] * delta;
      } else if(qii_m_qjj < 0) {
        S_sym[i + n*j] = p_0[j] * (expm1(qii_m_qjj*delta)/qii_m_qjj);
      } else { // ensure use of neg sign, avoids overflow in expm1
        S_sym[i + n*j] = p_0[i] * (expm1(-qii_m_qjj*delta)/(-qii_m_qjj));
      }
    }
  }
  
  /* populate X=S-Id */
  for( j=0; j<n; j++){
    for( i=0; i<n; i++){
      if(i==j){ // diagonal
        X[i + n*j] = expm1(delta*dQ[i]); // implicitly substracting 1
      } else if(i<j){ // upper triangle
        X[i + n*j] = Q[i + n*j] * S_sym[i + n*j];
      } else{ // lower triangle
        X[i + n*j] = Q[i + n*j] * S_sym[j + n*i]; // transpose S_sym
      }
    }
  }
  free(dQ);
  free(p_0);
  free(S_sym);
}

/* Implicit squaring for S=Id+X by looping X <- 2X + X^2 */
void matsq_hp(int n, int k, double *X)
{
  int i, inSqr=n*n, ione = 1;
  static const char *transN = "N";
  static const double one = 1.0, zero = 0.0, two = 2.0;
  size_t nSqr = n * ((size_t) n);
  double* w = (double*)malloc(sizeof(double) * nSqr); // temp storage for dgemm
  
  Memcpy(w, X, nSqr); // copy w <- X
  for (i = 0; i < k; i++){
    F77_CALL(dgemm)(transN, transN, &n, &n, &n, &one,
             X, &n, X, &n, &two, w, &n FCONE FCONE); // w <- 2w + X^2 = 2X+X^2
    Memcpy(X, w, nSqr); // X <- w = 2X+X^2
  }
  free(w);
}

// /* Implicit squaring for S=Id+X by looping X <- 2X + X^2 */
// void matsq_hp(int n, int k, double *X)
// {
//   int i, inSqr=n*n, ione = 1;
//   static const char *transN = "N";
//   static const double one = 1.0, zero = 0.0, two = 2.0;
//   size_t nSqr = n * ((size_t) n);
//   double* w = (double*)malloc(sizeof(double) * nSqr); // temp storage for dgemm
//   
//   for (i = 0; i < k; i++){
//     F77_CALL(dgemm)(transN, transN, &n, &n, &n, &one,
//              X, &n, X, &n, &zero, w, &n FCONE FCONE); // w <- X^2
//     F77_CALL(daxpy)(&inSqr, &two, X, &ione, w, &ione); // w <- w + 2X = X^2 + 2X
//     Memcpy(X, w, nSqr); // X <- w 
//   }
//   free(w);
// }
// 
// /* Regular squaring: compute z := z^{(2^K)}, via looping: w<-z*z; z<-w; (w is temp storage) */
// void matsq(int n, int k, double *z)
// {
//   int i;
//   static const char *transa = "N";
//   static const double one = 1.0, zero = 0.0;
//   size_t nSqr = n * ((size_t) n);
//   double* w = (double*)malloc(sizeof(double) * nSqr); // temp storage for dgemm
//   
//   for (i = 0; i < k; i++){
//     /* w <- z*z  */
//     F77_CALL(dgemm)(transa, transa, &n, &n, &n, &one,
//              z, &n, z, &n, &zero, w, &n FCONE FCONE);
//     Memcpy(z, w, nSqr); /* z <- w  */
//   }
//   
//   free(w);
// }
