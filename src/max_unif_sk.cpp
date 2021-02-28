#include "ProbSKE.hpp"
#include "ProbSKE.h" // for the dense versions of the methods
#include <thread>

/* 
 * max_sk_unif_expm: point-wise max of skeletoid and uniformization
 */

// .Call() this from R
RcppExport SEXP R_max_sk_unif_expm(SEXP Q_de_, SEXP Q_sp_, SEXP t_pow_,SEXP r_,
                                   SEXP n_, SEXP K_u_, SEXP K_s_, SEXP un_usp_,
                                   SEXP sk_usp_) {
  // Rcout << "Entered R_max_sk_unif_expm" << std::endl;
  const MapSpMatd Q_sp(as<MapSpMatd>(Q_sp_));
  const double t_pow(as<double>(t_pow_));
  const double r(as<double>(r_));
  const size_t n(as<size_t>(n_));
  const int K_u(as<int>(K_u_));
  const int K_s(as<int>(K_s_));
  const bool un_usp(as<bool>(un_usp_));
  const bool sk_usp(as<bool>(sk_usp_));
  MatrixXd ans(n,n);
  // Rcout << "Calling max_sk_unif_expm" << std::endl;
  max_sk_unif_expm(ans,REAL(Q_de_),Q_sp,t_pow,r,n,K_u,K_s,un_usp,sk_usp);
  return(wrap(ans));
}

void max_sk_unif_expm(MatrixXd& ans, double *Q_de, const MapSpMatd& Q_sp,
                      const double t_pow, const double r, const size_t n,
                      const int K_u,const int K_s, const bool un_usp,
                      const bool sk_usp){
  const double delta_s = t_pow*pow(2.0,-K_s);
  if(un_usp && sk_usp){
    // Rcout << "Both use sparse methods" << std::endl;
    MatrixXd E_u(n,n);
    MatrixXd E_s(n,n);
    std::thread t_u(unif_expm_sparse, std::ref(E_u), Q_sp, t_pow, r, n, K_u);
    std::thread t_s(skeletoid_expm_sparse, std::ref(E_s), Q_sp, delta_s, n, K_s);
    t_u.join();
    t_s.join();
    ans = E_u.array().max(E_s.array()).matrix();
  } else if(un_usp && !sk_usp){
    // Rcout << "Uniformization sparse and Skeletoid dense" << std::endl;
    MatrixXd E_u(n,n);
    double* E_s = (double*)malloc(sizeof(double) * n * n);
    std::thread t_u(unif_expm_sparse, std::ref(E_u), Q_sp, t_pow, r, n, K_u);
    std::thread t_s(skeletoid_expm_dense, E_s, Q_de, delta_s, n, K_s);
    t_u.join();
    t_s.join();
    ans = E_u.array().max(MapMatd(E_s,n,n).array()).matrix();
    free(E_s);
  } else if(!un_usp && sk_usp){
    // Rcout << "Uniformization dense and Skeletoid sparse" << std::endl;
    double* E_u = (double*)malloc(sizeof(double) * n * n);
    MatrixXd E_s(n,n);
    std::thread t_u(unif_expm_dense, E_u, Q_de, t_pow, r, n, K_u);
    std::thread t_s(skeletoid_expm_sparse, std::ref(E_s), Q_sp, delta_s, n, K_s);
    t_u.join();
    t_s.join();
    ans = E_s.array().max(MapMatd(E_u,n,n).array()).matrix();
    free(E_u);
  } else{
    // Rcout << "Both use dense methods" << std::endl;
    double* E_u = (double*)malloc(sizeof(double) * n * n);
    std::thread t_u(unif_expm_dense, E_u, Q_de, t_pow, r, n, K_u);
    double* E_s = (double*)malloc(sizeof(double) * n * n);
    std::thread t_s(skeletoid_expm_dense, E_s, Q_de, delta_s, n, K_s);
    t_u.join();
    t_s.join();
    ans = MapMatd(E_s,n,n).array().max(MapMatd(E_u,n,n).array()).matrix();
    free(E_u);
    free(E_s);
  }
}

/* 
 * max_sk_unif_vtexpm: point-wise max of skeletoid and uniformization
 */

// .Call() this from R
RcppExport SEXP R_max_sk_unif_vtexpm(SEXP Q_de_, SEXP Q_sp_, SEXP t_pow_,SEXP r_,
                                     SEXP nr_, SEXP nc_, SEXP K_u_, SEXP K_s_,
                                     SEXP k1_, SEXP k2_, SEXP v_de_, SEXP v_sp_,
                                     SEXP un_usp_, SEXP sk_usp_, SEXP full_sp_) {
  const MapSpMatd Q_sp(as<MapSpMatd>(Q_sp_));
  const double t_pow(as<double>(t_pow_));
  const double r(as<double>(r_));
  const size_t nr(as<size_t>(nr_));
  const size_t nc(as<size_t>(nc_));
  const int K_u(as<int>(K_u_));
  const int K_s(as<int>(K_s_));
  const int k1(as<int>(k1_));
  const int k2(as<int>(k2_));
  const MapSpMatd v_sp(as<MapSpMatd>(v_sp_));
  const bool un_usp(as<bool>(un_usp_));
  const bool sk_usp(as<bool>(sk_usp_));
  const bool full_sp(as<bool>(full_sp_));
  MatrixXd ans(nr,nc);
  max_sk_unif_vtexpm(ans,REAL(Q_de_),Q_sp,t_pow,r,nr,nc,K_u,K_s,k1,k2,
                     REAL(v_de_),v_sp,un_usp,sk_usp,full_sp);
  return(wrap(ans));
}

void max_sk_unif_vtexpm(MatrixXd& ans, double *Q_de, const MapSpMatd& Q_sp,
                        const double t_pow, const double r, const size_t nr,
                        const size_t nc, const int K_u,const int K_s, const int k1,
                        const int k2, double *v_de, const MapSpMatd& v_sp,
                        const bool un_usp,const bool sk_usp, const bool full_sp){
  const double delta_s = t_pow*pow(2.0,-K_s);
  if(un_usp && sk_usp){
    // Rcout << "Both use sparse methods" << std::endl;
    MatrixXd E_u(nr,nc);
    MatrixXd E_s(nr,nc);
    std::thread t_u(unif_vtexpm_sparse, std::ref(E_u), Q_sp, t_pow, r, v_sp,
                    nc, K_u);
    std::thread t_s(skeletoid_vtexpm_sparse, std::ref(E_s), Q_sp, delta_s, 
                    v_sp, nc, K_s, k1, k2, full_sp);
    t_u.join();
    t_s.join();
    ans = E_u.array().max(E_s.array()).matrix();
  } else if(un_usp && !sk_usp){
    // Rcout << "Uniformization sparse and Skeletoid dense" << std::endl;
    MatrixXd E_u(nr,nc);
    double* E_s = (double*)malloc(sizeof(double) * nr * nc);
    std::thread t_u(unif_vtexpm_sparse, std::ref(E_u), Q_sp, t_pow, r, v_sp,
                    nc, K_u);
    std::thread t_s(skeletoid_vtexpm_dense, E_s, Q_de, delta_s, v_de, nr, nc,
                    k1, k2);
    t_u.join();
    t_s.join();
    ans = E_u.array().max(MapMatd(E_s,nr,nc).array()).matrix();
    free(E_s);
  } else if(!un_usp && sk_usp){
    // Rcout << "Uniformization dense and Skeletoid sparse" << std::endl;
    double* E_u = (double*)malloc(sizeof(double) * nr * nc);
    MatrixXd E_s(nr,nc);
    std::thread t_u(unif_vtexpm_dense, E_u, Q_de, t_pow, r, nr, nc, K_u, v_de);
    std::thread t_s(skeletoid_vtexpm_sparse, std::ref(E_s), Q_sp, delta_s, 
                    v_sp, nc, K_s, k1, k2, full_sp);
    t_u.join();
    t_s.join();
    ans = E_s.array().max(MapMatd(E_u,nr,nc).array()).matrix();
    free(E_u);
  } else{
    // Rcout << "Both use dense methods" << std::endl;
    double* E_u = (double*)malloc(sizeof(double) * nr * nc);
    std::thread t_u(unif_vtexpm_dense, E_u, Q_de, t_pow, r, nr, nc, K_u, v_de);
    double* E_s = (double*)malloc(sizeof(double) * nr * nc);
    std::thread t_s(skeletoid_vtexpm_dense, E_s, Q_de, delta_s, v_de, nr, nc, k1, k2);
    t_u.join();
    t_s.join();
    ans = MapMatd(E_s,nr,nc).array().max(MapMatd(E_u,nr,nc).array()).matrix();
    free(E_u);
    free(E_s);
  }
}
