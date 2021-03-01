#include "pske.h"

/* 
 * unif_expm_dense: solution of the Kolmogorov's eq. via uniformization
 */

/* This is not called directly, see uniformization.cpp */
void unif_expm_dense_c(double *E, double *Q, double r, double *a, int n, int K){
  int i, ione = 1, inSqr = n*n;
  double inv_r = 1.0/r;
  set_diagonal(E,n,a[0]); // set E[i,i]=exp(-r) and 0 o.w.
  if(K==0) return;
  
  // define and initialize P = diag(n) + Q/r
  double* P = (double*)malloc(sizeof(double) * inSqr);
  uniformize_ratemat_dense(P,Q,inv_r,n,inSqr,ione);
  
  // define and initialize P_pow = P^1 = P
  double* P_pow = (double*)malloc(sizeof(double) * inSqr);
  Memcpy(P_pow, P, inSqr); // P_pow <- P
  uniformization_dense_loop(P, P_pow, a, n, n, K, E);
  free(P);
  free(P_pow);
}

/* runs the inner loop associated with the series expression */
void uniformization_dense_loop(double *P, double *P_pow, double *a, int nr,
                               int nc, int K, double *E){
  int i, ione = 1, inr_t_nc = nr*nc;
  static const char *transN = "N";
  static const double one = 1.0, zero = 0.0;
  double tmp;
  double* w = (double*)malloc(sizeof(double) * inr_t_nc); // temp storage for dgemm
  
  for (i = 1; i <= K; i++){
    tmp=a[i];
    F77_CALL(daxpy)(&inr_t_nc, &tmp, P_pow, &ione, E, &ione); // E <- a[i]*P_pow + E
    F77_CALL(dgemm)(transN, transN, &nr, &nc, &nc, &one, // w <- P_pow %*% P
             P_pow, &nr, P, &nc, &zero, w, &nr FCONE FCONE);
    Memcpy(P_pow, w, inr_t_nc); // P_pow <- w
  }
  free(w);
}

/* 
 * unif_vtexpm_dense: solution of the evolution eqs. via uniformization
 */

/* do the actual work */
void unif_vtexpm_dense_c(double *E, double *Q, double r, double *a, int nr,
                         int nc, int K, double *v){
  int i, ione = 1, inSqr = nc*nc, nr_t_nc = nr*nc;
  static const char *transN = "N";
  static const double one = 1.0, zero = 0.0;
  double inv_r = 1.0/r, tmp=a[0];
  
  // initialize E = a[0]*v
  Memcpy(E, v, nr_t_nc); // E <- v
  F77_CALL(dscal)(&nr_t_nc, &tmp, E, &ione); // E <- a[0]*E = a[0]*v
  if(K==0) return;
  
  // define and initialize P = diag(n) + Q/r
  double* P = (double*)malloc(sizeof(double) * inSqr);
  uniformize_ratemat_dense(P,Q,inv_r,nc,inSqr,ione);
  
  // define and initialize vP_pow = v*P^1 = vP
  double* vP_pow = (double*)malloc(sizeof(double) * nr_t_nc);
  F77_CALL(dgemm)(transN, transN, &nr, &nc, &nc, &one,
           v, &nr, P, &nc, &zero, vP_pow, &nr FCONE FCONE);
  
  // do inner loop
  uniformization_dense_loop(P, vP_pow, a, nr, nc, K, E);
  free(P);
  free(vP_pow);
}

// /* runs the inner loop associated with the series expression */
// void unif_vtexpm_dense_loop(double *P, double *vP_pow,
//                             double *a, int nr, int nc, int K, double *E)
// {
//   int i, ione = 1;
//   static const char *transT = "T";
//   static const double one = 1.0, zero = 0.0;
//   double tmp;
//   double* w = (double*)malloc(sizeof(double) * n); // temp storage for dgemv
//   
//   for (i = 1; i <= K; i++){
//     tmp = a[i];
//     F77_CALL(daxpy)(&n, &tmp, vP_pow, &ione, E, &ione); // E <- E + a[i]*vP_pow
//     F77_CALL(dgemv)(transT, &n, &n, &one, P, // w = vP_pow^T%*%P = (P^T%*%vP_pow)^T
//              &n, vP_pow, &ione, &zero, w, &ione FCONE);
//     Memcpy(vP_pow, w, n); /* P_pow <- w */
//   }
//   free(w);
// }

/* 
 * utils
 */

// write a diagonal matrix on the n*n array A with constant value a
void set_diagonal(double *A, int n, double a){
  for(int j=0;j<n;j++)
    for(int i=0;i<n;i++)
      A[i+n*j] = ((i==j)?(a):(0.0));
}

/* define and initialize P = diag(n) + Q/r
 * 1) create P = diag(n), then 2) use daxpy with P<-P+(1/r)*Q */ 
void uniformize_ratemat_dense(double *P, double *Q, double inv_r, int n, 
                         int inSqr, int ione){
  set_diagonal(P,n,1.0); // set P = identity(n)
  F77_CALL(daxpy)(&inSqr, &inv_r, Q, &ione, P, &ione); // P<-P+(1/r)*Q
}

// Find minimum along the diagonal of a matrix (O(n))
// For use when R's efficient min() (O(1)) is not available
double min_diag(double *A, int n){
  int p;
  double min = A[0];
  for(int i=1;i<n;i++){
    p=i+n*i;
    if(A[p]<min)
      min=A[p];
  }
  return(min);
}