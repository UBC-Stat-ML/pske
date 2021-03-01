##############################################################################
# full expm solver
##############################################################################

#' Uniformization
#'
#' @param Q rate matrix.
#' @param t_pow the t in exp(tQ)
#' @param eps desired error tolerance
#' @param K precision level. If provided, eps is ignored.
#' @returns exp(tQ) up to tolerance eps in op-1 norm.
#' @export
unif_expm = function(Q, K, t_pow=1, eps, sparse){
  r=get_max_rate(Q)
  if(missing(K)) K = unif_auto_tune(Q=Q,t_pow=t_pow,eps=eps,r=r)$K
  if(missing(sparse)) sparse = unif_expm_use_sparse(Q)
  if(sparse){
    return(unif_expm_sparse(Q, K, t_pow, r))
  } else{
    return(unif_expm_dense(Q, K, t_pow, r))
  }
}

###############################################################################
# helpers and wrappers to compiled code
###############################################################################

unif_expm_use_sparse = function(Q){
  (is(Q,"sparseMatrix") && Matrix::nnzero(Q) <= 0.075*(nrow(Q))^2)
}

# expm sparse case
unif_expm_sparse = function(Q, K, t_pow, r){
  .Call("R_unif_expm_sparse", PACKAGE = 'pske',
        as(Q,"dgCMatrix"), as.double(t_pow), as.double(r), nrow(Q),
        as.integer(K))
}

# expm dense case
unif_expm_dense = function(Q, K, t_pow, r){
  .Call("R_unif_expm_dense", PACKAGE = 'pske',
        as(Q,"matrix"), as.double(t_pow), as.double(r),
        nrow(Q), as.integer(K))
}