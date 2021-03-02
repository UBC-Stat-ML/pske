#include "pske.hpp"
#include "pske.h" // for the dense versions of uniformization

/* 
 * Uses Alg 1 in Sherlock (2020) for stability at high rates
 */

// .Call() this from R
RcppExport SEXP R_unif_expm_sparse(SEXP Q_, SEXP t_pow_, SEXP r_, SEXP n_,
                                   SEXP K_) {
  const MapSpMatd Q(as<MapSpMatd>(Q_));
  const double t_pow(as<double>(t_pow_));
  const double r(as<double>(r_));
  const size_t n(as<size_t>(n_));
  const int K(as<int>(K_));
  MatrixXd ans = MatrixXd::Identity(n,n);
  unif_sparse(ans,Q,t_pow,r,n,K);
  return(wrap(ans));
}

// .Call() this from R
RcppExport SEXP R_unif_vtexpm_sparse(SEXP Q_, SEXP t_pow_, SEXP r_,
                                     SEXP v_, SEXP n_, SEXP K_) {
  const MapSpMatd Q(as<MapSpMatd>(Q_));
  const MapSpMatd v(as<MapSpMatd>(v_));
  const double t_pow(as<double>(t_pow_));
  const double r(as<double>(r_));
  const size_t n(as<size_t>(n_));
  const int K(as<int>(K_));
  MatrixXd ans(v);
  unif_sparse(ans,Q,t_pow,r,n,K);
  return(wrap(ans));
}

void unif_sparse(MatrixXd& E, const MapSpMatd& Q, const double t_pow,
                      const double r, const size_t n, const int K){
  double lambda = r*t_pow; // get rate of the homogeneous Poisson process
  if(K == 0){
    E*=exp(-lambda);
    return;
  }
  
  // get M := rtP = Qt + rtI = Qt + lambdaI
  SpMatd M(t_pow*Q);
  M.diagonal().array() += lambda;
  
  // series loop
  MatrixXd M_pow(E); // initialize M_pow=E=Id as dense since it store powers of P
  uniformization_loop_sparse(E, M_pow, M, lambda, K);
}

// unif loop: uses Alg 1 in Sherlock (2020) for stability at high rates
// note: templating is very slow because an extra copy is made. See
// https://stackoverflow.com/a/57427407/5443023
// template <typename Derived>
void uniformization_loop_sparse(MatrixXd& E, MatrixXd& S,const SpMatd& M,
                                // const Eigen::EigenBase<Derived>& M_,
                                const double lambda, const int K){
  // Derived const& M = M_.derived(); // retrieve actual matrix (sparse or dense)
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

/* 
 * unif_dense
 */

// .Call() this from R
RcppExport SEXP R_unif_expm_dense(SEXP Q_, SEXP t_pow_, SEXP r_, SEXP n_,
                                  SEXP K_) {
  const MapMatd Q(as<MapMatd>(Q_));
  const double t_pow(as<double>(t_pow_));
  const double r(as<double>(r_));
  const size_t n(as<size_t>(n_));
  const int K(as<int>(K_));
  MatrixXd ans = MatrixXd::Identity(n,n);
  unif_dense(ans,Q,t_pow,r,n,K);
  return(wrap(ans));
}

// .Call() this from R
RcppExport SEXP R_unif_vtexpm_dense(SEXP Q_, SEXP t_pow_, SEXP r_,
                                    SEXP v_, SEXP n_, SEXP K_) {
  const MapMatd Q(as<MapMatd>(Q_));
  const MapMatd v(as<MapMatd>(v_));
  const double t_pow(as<double>(t_pow_));
  const double r(as<double>(r_));
  const size_t n(as<size_t>(n_));
  const int K(as<int>(K_));
  MatrixXd ans(v);
  unif_dense(ans,Q,t_pow,r,n,K);
  return(wrap(ans));
}

void unif_dense(MatrixXd& E, const MapMatd& Q, const double t_pow,
                const double r, const size_t n, const int K){
  double lambda = r*t_pow; // get rate of the homogeneous Poisson process
  if(K == 0){
    E*=exp(-lambda);
    return;
  }
  
  // get M := rtP = Qt + rtI = Qt + lambdaI
  MatrixXd M(t_pow*Q);
  M.diagonal().array() += lambda;
  
  // series loop
  MatrixXd M_pow(E); // initialize M_pow=E=Id as dense since it store powers of P
  uniformization_loop_dense(E, M_pow, M, lambda, K);
}

// unif loop: uses Alg 1 in Sherlock (2020) for stability at high rates
// template <typename Derived>
void uniformization_loop_dense(MatrixXd& E, MatrixXd& S,const MatrixXd& M,
                            // const Eigen::EigenBase<Derived>& M_,
                               const double lambda, const int K){
  // Derived const& M = M_.derived(); // retrieve actual matrix (sparse or dense)
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