#include "ProbSKE.h"

/* .Call() this from R */
SEXP R_skeletoid_vtexpm_dense(SEXP Q_, SEXP delta_, SEXP v_,
                              SEXP nr_, SEXP nc_, SEXP k1_, SEXP k2_) {
  R_len_t nr=asInteger(nr_),nc=asInteger(nc_),k1=asInteger(k1_),k2=asInteger(k2_);
  double delta = asReal(delta_);
  int nprotect = 0;
  SEXP A_sexp; // declare answer
  PROTECT(A_sexp = allocMatrix(REALSXP, nr, nc)); nprotect++; //allocate+PROTECT
  skeletoid_vtexpm_dense(REAL(A_sexp), REAL(Q_), delta, REAL(v_),
                         nr, nc, k1, k2);
  UNPROTECT(nprotect);
  return(A_sexp);
}

/* calls skeletoid_expm_dense for k1 squarings, then does 2^k2 vec*mat ops */
void skeletoid_vtexpm_dense(double *A, double *Q, double delta, 
                            double *v, size_t nr, size_t nc, int k1, int k2){
  double* S = (double*)malloc(sizeof(double) * nc * nc); // allocate space for S_delta
  skeletoid_expm_dense(S,Q,delta,nc,k1); // populate S=S_delta^(2^k1)
  int nvmops = 1 << k2; // 2^(k2) using bitshift, faster & +stable than (int)pow(2,k2)
  mat_mat_prod_loop(A, S, v, nr, nc, nvmops, true); // do init A=v*S; then loop: A<-A*S;
  free(S);
}

/* init: z=v*x; loop: z<-z*x; with v (nr times nc) and x (nc times nc) */
void mat_mat_prod_loop(double *z, double *x, double *v, int nr, int nc, int k,
                       bool do_init)
{
  int i, ione = 1;
  size_t nr_t_nc = nr * nc;
  static const char *transN = "N";
  static const double one = 1.0, zero = 0.0;
  double* w = (double*)malloc(sizeof(double) * nr_t_nc); // temp storage for dgemm
  if(do_init){
    // initialize z<-v*x
    F77_CALL(dgemm)(transN, transN, &nr, &nc, &nc, &one,
             v, &nr, x, &nc, &zero, z, &nr FCONE FCONE);
    k-=1; // subtract this initialization round
  }
  for (i = 0; i < k; i++){
    // w <- z*x
    F77_CALL(dgemm)(transN, transN, &nr, &nc, &nc, &one,
             z, &nr, x, &nc, &zero, w, &nr FCONE FCONE);
    Memcpy(z, w, nr_t_nc); // z <- w
  }
  free(w);
}

// /* init: z=v*x; loop: z<-z*x; */
// void vec_mat_prod_loop(double *z, double *x, double *v, int n, int k, bool do_init)
// {
//   int i, ione = 1;
//   static const char *transT = "T";
//   static const double one = 1.0, zero = 0.0;
//   
//   if(do_init){
//     // initialize z= v^T %*% x = (x^T %*% v)^T
//     F77_CALL(dgemv)(transT, &n, &n, &one, x,
//              &n, v, &ione, &zero, z, &ione FCONE);
//     k-=1; // subtract this initialization round
//   }
//   
//   double* w = (double*)malloc(sizeof(double) * n); // temp storage for dgemv
//   for (i = 0; i < k; i++){
//     /* w = z^T %*% x = (x^T %*% z)^T */
//     F77_CALL(dgemv)(transT, &n, &n, &one, x,
//              &n, z, &ione, &zero, w, &ione FCONE);
//     Memcpy(z, w, n); /* z <- w  */
//   }
//   free(w);
// }