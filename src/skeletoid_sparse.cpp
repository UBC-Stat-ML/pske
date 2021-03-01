#include "pske.hpp"
#include "pske.h" // for matsq

/* 
 * expm functions
 */

/* .Call() this from R */
RcppExport SEXP R_skeletoid_expm_sparse(SEXP Q_, SEXP delta_, SEXP n_sq_) {
  const SpMatd Q(as<MapSpMatd>(Q_)); // input: sparse matrix
  const double delta(as<double>(delta_)); // delta=t_pow*2^K
  const int n_sq(as<int>(n_sq_)); // squaring steps
  const size_t n = Q.rows();
  MatrixXd ans(n,n);
  skeletoid_expm_sparse(ans,Q,delta,n,n_sq);
  return(wrap(ans));
}

/* build S_delta sparse end then compute S_delta^2^n_sq using dense ops */
void skeletoid_expm_sparse(MatrixXd& E, const SpMatd& Q, const double delta,
                           const size_t n, const int n_sq) {
  SpMatd X_delta(n,n); // define X_delta = S_delta - Id
  build_S_delta_sparse(X_delta,Q,delta,n,true); // populate X_delta
  E = MatrixXd(X_delta); // cast to dense
  matsq_hp(n,n_sq,E.data()); // squaring loop using direct BLAS calls
  E.diagonal().array() += 1.0; // Recover S^(2^n_sq) = Id + X^(2^n_sq)
}

/* Takes an allocated S_delta sparse matrix and fills it */
void build_S_delta_sparse(SpMatd& S_delta, const SpMatd& Q, const double delta,
                          const size_t n, const bool hp) {
  const size_t nnz = Q.nonZeros();
  ArrayXd d(Q.diagonal()); // extract diag(Q)
  ArrayXd p0((d*delta).exp()); // exp(delta*diag(Q))
  
  // compute elements of S_delta as triplets
  size_t i,j,k;
  double dd,val;
  std::vector<Tripd> tripletList;
  tripletList.reserve(nnz);
  for (k=0; k<Q.outerSize(); ++k)
    for (SparseMatrix<double>::InnerIterator it(Q,k); it; ++it){
      i=it.row(); j=it.col();
      if(i != j){ // need to do diagonal separately bc we could have qii=0, which would be skipped here
        dd = d(i) - d(j);
        if(abs(dd)<100*DBL_EPSILON){
          val=p0(i)*delta;
        } else if(dd<0){
          val=p0(j)*expm1(delta*dd)/dd;
        } else{
          val=p0(i)*expm1(-delta*dd)/(-dd);
        }
        tripletList.push_back(Tripd(i,j,it.value()*val)); // fill with q_{ij} * val
      }
    }
  // fill diagonal
  if(hp)
    for (i=0; i<n; i++){
      val = expm1(d(i)*delta); // TODO: do vectorized expm1 when adopted in Eigen stable
      tripletList.push_back(Tripd(i,i,val));
    }
  else
    for (i=0; i<n; i++)
      tripletList.push_back(Tripd(i,i,p0(i)));
  S_delta.setFromTriplets(tripletList.begin(), tripletList.end()); // populate sparse matrix
}

/* 
 * vtexpm functions
 */

/* .Call() this from R*/
RcppExport SEXP R_skeletoid_vtexpm_sparse(SEXP Q_, SEXP delta_, SEXP v_,
                                          SEXP K_, SEXP k1_, SEXP k2_,
                                          SEXP full_sp_) {
  const SpMatd Q(as<MapSpMatd>(Q_)); // input: sparse matrix
  const SpMatd v(as<MapSpMatd>(v_)); // input: sparse matrix (rowvector)
  const double delta(as<double>(delta_)); // delta=t_pow*2^K
  const int K(as<int>(K_)); // number of vec*mat ops
  const int k1(as<int>(k1_)); // number of squarings
  const int k2(as<int>(k2_)); // number of vec*mat ops
  const bool full_sp(as<bool>(full_sp_)); // all sparse ops or just building S_delta
  const size_t nc = Q.cols();
  MatrixXd ans(v.rows(),nc);
  skeletoid_vtexpm_sparse(ans,Q,delta,v,nc,K,k1,k2,full_sp);
  return(wrap(ans));
}

/* Builds S_delta sparse and then does either full sparse or mixed computations
 * In both cases we use A dense since nonzeros propagate fast in A
 */
void skeletoid_vtexpm_sparse(MatrixXd& A, const SpMatd& Q,
                             const double delta, const SpMatd& v,
                             const size_t nc, const int K, const int k1,
                             const int k2, const bool full_sp){
  SpMatd S_delta(nc,nc);
  const int nvmops = (1 << k2) - 1; // 2^(k2) - 1 = number of mat*mat ops (not counting init)
  // const int nvmops = (int)(pow(2.0,k2))-1; // number of mat*mat ops (not counting init)
  if(full_sp){
    build_S_delta_sparse(S_delta,Q,delta,nc,false); // no need for hp since only mat*mat ops
    skeletoid_vtexpm_sparse_full(A,S_delta,v,K,nvmops);
  } else{
    build_S_delta_sparse(S_delta,Q,delta,nc,true); // need hp because there might be squaring
    skeletoid_vtexpm_sparse_mixed(A,S_delta,v,nc,k1,k2,nvmops);
  }
}

/* does only dense*sparse mat ops */
void skeletoid_vtexpm_sparse_full(MatrixXd& A, SpMatd& S_delta,
                                  const SpMatd& v, const int K, const int nvmops){
  A = v*S_delta; // initialize A (1st mat*mat op, not counted below)
  for(int i=0; i<nvmops; i++)
    A*=S_delta;
}

/* Casts to dense and does almost all the work through calls to BLAS */
void skeletoid_vtexpm_sparse_mixed(MatrixXd& A, SpMatd& X_delta,
                                   const SpMatd& v, const size_t nc,
                                   const int k1, const int k2, const int nvmops){
  MatrixXd E = MatrixXd(X_delta); // cast X_delta to dense
  if(k1>0) 
    matsq_hp(nc,k1,E.data()); // implicit squaring loop for S via X<-2X+X^2
  E.diagonal().array() += 1.0; // Recover S^(2^n_sq) = Id + X^(2^n_sq)
  A = v*E; // initialize A (1st mat*mat op, not counted below)
  mat_mat_prod_loop(A.data(), E.data(), A.data(), A.rows(), nc, nvmops, false); // no re-init, loop: A<-A*E;
}
