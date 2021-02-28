###############################################################################
# skeletoid: full matrix exponential solver
###############################################################################

#' Skeletoid
#'
#' @param Q rate matrix.
#' @param t_pow the t in exp(tQ)
#' @param eps desired error tolerance
#' @param K precision level. If provided, eps is ignored.
#' @returns exp(tQ) up to tolerance eps in op-1 norm.
#' @export
skeletoid_expm = function(Q, K, t_pow=1, eps, n_sq = K, sparse){
  if(missing(K)) K = skeletoid_auto_tune(Q = Q,t_pow = t_pow,eps = eps)$K
  if(missing(sparse)) sparse = skeletoid_expm_use_sparse(Q)
  delta = t_pow*(2^(-K)) # size of increment
  if(sparse)
    return(skeletoid_expm_sparse(Q, delta, n_sq))
  else
    return(skeletoid_expm_dense(as(Q,"matrix"),delta,n_sq))
}

###############################################################################
# skeletoid: helpers and wrappers to compiled code
###############################################################################

# should we take advantage of sparse computations?
skeletoid_expm_use_sparse = function(Q){
  (is(Q,"sparseMatrix") && Matrix::nnzero(Q)<=0.2*(nrow(Q))^2) # from: test_switch_sparse_dense.R
}

# skeletoid: dense case
skeletoid_expm_dense = function(Q, delta, n_sq){
  .Call("R_skeletoid_expm_dense", PACKAGE = 'ProbSKE',
        as.double(Q), as.double(delta), nrow(Q), as.integer(n_sq))
}

# skeletoid: sparse case
skeletoid_expm_sparse = function(Q, delta, n_sq){
  .Call('R_skeletoid_expm_sparse', PACKAGE = 'ProbSKE',
        as(Q,"dgCMatrix"), as.double(delta), as.integer(n_sq))
}
