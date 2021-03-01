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
RcppExport SEXP R_unif_expm_sparse(SEXP Q_, SEXP t_pow_, SEXP r_, SEXP n_,
                                   SEXP K_);
void unif_expm_sparse(MatrixXd& E, const MapSpMatd& Q, const double t_pow,
                      const double r, const size_t n, const int K);
void get_Poisson_pmf(double *a, const double lambda, const int K);
void uniformize_ratemat_sparse(SpMatd& P,const MapSpMatd& Q,const double r);
void uniformization_loop_sparse(MatrixXd& E, MatrixXd& S, const SpMatd& P,
                                double *a, const int K);
RcppExport SEXP R_unif_vtexpm_sparse(SEXP Q_, SEXP t_pow_, SEXP r_,
                                     SEXP v_, SEXP n_, SEXP K_);
void unif_vtexpm_sparse(MatrixXd& E, const MapSpMatd& Q, const double t_pow,
                        const double r, const MapSpMatd& v, const size_t n,
                        const int K);
RcppExport SEXP R_unif_expm_dense(SEXP Q_, SEXP t_pow_, SEXP r_, SEXP n_, 
                                  SEXP K_);
void unif_expm_dense(double *E, double *Q, const double t_pow,
                     const double r, const size_t n, const int K);
RcppExport SEXP R_unif_vtexpm_dense(SEXP Q_, SEXP t_pow_, SEXP r_, SEXP nr_,
                                    SEXP nc_, SEXP K_, SEXP v_);
void unif_vtexpm_dense(double *E, double *Q, const double t_pow,
                       const double r, const size_t nr, const size_t nc,
                       const int K, double *v);

// functions in max_sk_unif.cpp
RcppExport SEXP R_max_sk_unif_expm(SEXP Q_de_, SEXP Q_sp_, SEXP t_pow_,SEXP r_,
                                   SEXP n_, SEXP K_u_, SEXP K_s_, SEXP un_usp_,
                                   SEXP sk_usp_);
void max_sk_unif_expm(MatrixXd& ans, double *Q_de, const MapSpMatd& Q_sp,
                        const double t_pow, const double r, const size_t n,
                        const int K_u,const int K_s, const bool un_usp,
                        const bool sk_usp);
RcppExport SEXP R_max_sk_unif_vtexpm(SEXP Q_de_, SEXP Q_sp_, SEXP t_pow_,SEXP r_,
                                     SEXP nr_, SEXP nc_, SEXP K_u_, SEXP K_s_,
                                     SEXP k1_, SEXP k2_, SEXP v_de_, SEXP v_sp_,
                                     SEXP un_usp_, SEXP sk_usp_, SEXP full_sp_);
void max_sk_unif_vtexpm(MatrixXd& ans, double *Q_de, const MapSpMatd& Q_sp,
                        const double t_pow, const double r, const size_t nr,
                        const size_t nc, const int K_u,const int K_s, const int k1,
                        const int k2, double *v_de, const MapSpMatd& v_sp,
                        const bool un_usp,const bool sk_usp, const bool full_sp);

#endif
