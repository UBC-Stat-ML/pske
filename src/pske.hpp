#ifndef PSKE_HPP
#define PSKE_HPP
#include <RcppEigen.h>

using Rcpp::as;
using Rcpp::wrap;
using Rcpp::Rcout;
using Eigen::Map;
using Eigen::Matrix;
using Eigen::SparseMatrix;
using Eigen::VectorXi;
using Eigen::ArrayXd;
using Eigen::MatrixXd;
typedef Map<MatrixXd> MapMatd;
typedef Map<SparseMatrix<double>> MapSpMatd;
typedef SparseMatrix<double> SpMatd;
typedef Eigen::Triplet<double> Tripd;

/* Note on the "RcppExport" tag: it is essentially == extern "C"
 * And the idea of using extern "C" here is to be able to call C++ code from C
 * code, which the R interpreter is based on. More info:
 * https://stackoverflow.com/a/30526795/5443023
 */

// functions in skeletoid_sparse.cpp
RcppExport SEXP R_skeletoid_expm_sparse(SEXP Q_, SEXP delta_, SEXP n_sq_);
void skeletoid_expm_sparse(MatrixXd& E, const SpMatd& Q, const double delta,
                           const size_t n, const int n_sq);
void build_S_delta_sparse(SpMatd& S_delta, const SpMatd& Q, const double delta,
                          const size_t n, const bool hp);
RcppExport SEXP R_skeletoid_vtexpm_sparse(SEXP Q_, SEXP delta_, SEXP v_,
                                          SEXP K_, SEXP k1_, SEXP k2_,
                                          SEXP full_sp_);
void skeletoid_vtexpm_sparse(MatrixXd& A, const SpMatd& Q,
                             const double delta, const SpMatd& v,
                             const size_t nc, const int K, const int k1,
                             const int k2, const bool full_sp);
void skeletoid_vtexpm_sparse_full(MatrixXd& A, SpMatd& S_delta,
                                  const SpMatd& v, const int K, const int nvmops);
void skeletoid_vtexpm_sparse_mixed(MatrixXd& A, SpMatd& X_delta,
                                   const SpMatd& v, const size_t nc,
                                   const int k1, const int k2, const int nvmops);

// functions in uniformization.cpp
// sparse case
RcppExport SEXP R_unif_expm_sparse(SEXP Q_, SEXP t_pow_, SEXP r_, SEXP n_,
                                   SEXP K_);
RcppExport SEXP R_unif_vtexpm_sparse(SEXP Q_, SEXP t_pow_, SEXP r_,
                                     SEXP v_, SEXP n_, SEXP K_);
void unif_sparse(MatrixXd& E, const MapSpMatd& Q, const double t_pow,
                 const double r, const size_t n, const int K);
// dense case
RcppExport SEXP R_unif_expm_dense(SEXP Q_, SEXP t_pow_, SEXP r_, SEXP n_,
                                  SEXP K_);
RcppExport SEXP R_unif_vtexpm_dense(SEXP Q_, SEXP t_pow_, SEXP r_,
                                    SEXP v_, SEXP n_, SEXP K_);
void unif_dense(MatrixXd& E, const MapMatd& Q, const double t_pow,
                const double r, const size_t n, const int K);

// templated functions for uniformization.cpp
// unif loop: uses Alg 1 in Sherlock (2020) for stability at high rates
template <typename Derived>
void uniformization_loop(MatrixXd& E, MatrixXd& S, 
                         const Eigen::EigenBase<Derived>& M_,
                         const double lambda, const int K){
  Derived const& M = M_.derived(); // retrieve actual matrix (sparse or dense)
  double b=S.lpNorm<1>(), c=0.0;
  const double B = 1e100;
  if(b>B){ // need to renormalize
    S = S/b; E = E/b; c += log(b); b=1.0;
  }
  for(int k=1; k<=K; k++){ // we already did i=0 when initializing E
    S = S*M/k; // update the product (dense*sparse)
    b = b*lambda/k; // update norm estimate
    E = E + S; // accumulate
    if(b>B){ // need to renormalize
      S = S/b; E = E/b; c += log(b); b=1.0;
    }
  }
  E = exp(c-lambda)*E; // undo renormalizations
  return;
}


#endif
