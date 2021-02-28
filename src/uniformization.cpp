#include "ProbSKE.hpp"
#include "ProbSKE.h" // for the dense versions of uniformization
#include <boost/math/distributions/poisson.hpp>

/* 
 * unif_expm_sparse: solution of the Kolmogorov's eq. via uniformization
 */

// .Call() this from R
RcppExport SEXP R_unif_expm_sparse(SEXP Q_, SEXP t_pow_, SEXP r_, SEXP n_,
                                   SEXP K_) {
  const MapSpMatd Q(as<MapSpMatd>(Q_));
  const double t_pow(as<double>(t_pow_));
  const double r(as<double>(r_));
  const size_t n(as<size_t>(n_));
  const int K(as<int>(K_));
  MatrixXd ans(n,n);
  unif_expm_sparse(ans,Q,t_pow,r,n,K);
  return(wrap(ans));
}

void unif_expm_sparse(MatrixXd& E, const MapSpMatd& Q, const double t_pow,
                      const double r, const size_t n, const int K){
  // compute Poisson pmf's for i=0...K
  // double r = -(Q.coeffs().minCoeff());
  double lambda = r*t_pow; // get rate of the homogeneous Poisson process
  double* a = (double*)malloc(sizeof(double) * (K+1)); // declare array of pmfs
  get_Poisson_pmf(a,lambda,K); // fill a
  
  // initialize accumulator
  E.setIdentity();
  E*=a[0]; // now E = exp(-r)*diag(n), which is case K=0
  if(K == 0) return; // so we're done if K=0
  
  // uniformize Q into P
  SpMatd P(n,n); // declare transition probability matrix
  uniformize_ratemat_sparse(P,Q,r);
  
  // series loop
  MatrixXd P_pow(P); // initialize P_pow=P^1 as dense since it store powers of P
  uniformization_loop_sparse(E, P_pow, P, a, K);
  free(a);
}

void get_Poisson_pmf(double *a, const double lambda, const int K){
  boost::math::poisson_distribution<double> pois(lambda); // instantiate Poisson dist
  for(int i=0;i<=K;i++)
    a[i]=pdf(pois,i);
}

void uniformize_ratemat_sparse(SpMatd& P,const MapSpMatd& Q,const double r){
  P.setIdentity(); // initialize P to identity
  P = P + Q/r; // P <- diag(n) + Q/r
}

void uniformization_loop_sparse(MatrixXd& E, MatrixXd& S, const SpMatd& P,
                                double *a, const int K){
  for(int i=1; i<=K; i++){ // we already did i=0 when initializing E
    E = E + a[i]*S;
    S = S*P; // dense*sparse product so we still take advantage of sparsity
  }
}

/* 
 * unif_vtexpm_sparse: solution of the evolution equation via uniformization
 */

// .Call() this from R
RcppExport SEXP R_unif_vtexpm_sparse(SEXP Q_, SEXP t_pow_, SEXP r_,
                                     SEXP v_, SEXP n_, SEXP K_) {
  const MapSpMatd Q(as<MapSpMatd>(Q_));
  const MapSpMatd v(as<MapSpMatd>(v_));
  const double t_pow(as<double>(t_pow_));
  const double r(as<double>(r_));
  const size_t n(as<size_t>(n_));
  const int K(as<int>(K_));
  MatrixXd ans(v.rows(),n);
  unif_vtexpm_sparse(ans,Q,t_pow,r,v,n,K);
  return(wrap(ans));
}

void unif_vtexpm_sparse(MatrixXd& E, const MapSpMatd& Q, const double t_pow,
                        const double r, const MapSpMatd& v, const size_t n,
                        const int K){
  // double r = -(Q.coeffs().minCoeff()); // get maximal rate
  double lambda = r*t_pow; // get rate of the homogeneous Poisson process
  double* a = (double*)malloc(sizeof(double) * (K+1)); // declare array of pmfs
  get_Poisson_pmf(a,lambda,K); // fill a
  
  // initialize accumulator
  E = a[0]*v; // now E = exp(-r)*v, which is case K=0
  if(K == 0) return; // so we can safely return if K=0
  
  // uniformize Q
  SpMatd P(n,n); // declare transition probability matrix
  uniformize_ratemat_sparse(P,Q,r); // fill P
  
  // series loop
  MatrixXd vP_pow(v*P); // initialize vP_pow=v*P^1 as dense since it store powers of P
  uniformization_loop_sparse(E, vP_pow, P, a, K);
  free(a);
}

/* 
 * unif_expm_dense: thin wrapper on the C code, provides Poisson pmfs
 */

// .Call() this from R 
RcppExport SEXP R_unif_expm_dense(SEXP Q_, SEXP t_pow_, SEXP r_, SEXP n_, 
                                  SEXP K_) {
  const double t_pow(as<double>(t_pow_));
  const double r(as<double>(r_));
  const size_t n(as<size_t>(n_));
  const int K(as<int>(K_));
  MatrixXd ans(n,n);
  unif_expm_dense(ans.data(), REAL(Q_), t_pow, r, n, K);
  return(wrap(ans));
}

void unif_expm_dense(double *E, double *Q, const double t_pow,
                     const double r, const size_t n, const int K){
  // compute Poisson pmf's for i=0...K
  double lambda = r*t_pow; // get rate of the homogeneous Poisson process
  double* a = (double*)malloc(sizeof(double) * (K+1)); // declare array of pmfs
  get_Poisson_pmf(a,lambda,K); // fill a
  unif_expm_dense_c(E,Q,r,a,n,K); // call C code
  free(a);
}

/* 
 * unif_vtexpm_dense: thin wrapper on the C code, provides Poisson pmfs
 */

// .Call() this from R 
RcppExport SEXP R_unif_vtexpm_dense(SEXP Q_, SEXP t_pow_, SEXP r_, SEXP nr_,
                                    SEXP nc_, SEXP K_, SEXP v_) {
  const double t_pow(as<double>(t_pow_));
  const double r(as<double>(r_));
  const size_t nr(as<size_t>(nr_));
  const size_t nc(as<size_t>(nc_));
  const int K(as<int>(K_));
  int nprotect = 0;
  MatrixXd ans(nr,nc);
  unif_vtexpm_dense(ans.data(), REAL(Q_), t_pow, r, nr, nc, K, REAL(v_));
  return(wrap(ans));
}

void unif_vtexpm_dense(double *E, double *Q, const double t_pow,
                       const double r, const size_t nr, const size_t nc,
                       const int K, double *v){
  // compute Poisson pmf's for i=0...K
  double lambda = r*t_pow; // get rate of the homogeneous Poisson process
  double* a = (double*)malloc(sizeof(double) * (K+1)); // declare array of pmfs
  get_Poisson_pmf(a,lambda,K); // fill a
  unif_vtexpm_dense_c(E,Q,r,a,nr,nc,K,v); // call C code
  free(a);
}
