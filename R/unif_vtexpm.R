##############################################################################
# uniformization evolution equation solver
##############################################################################

#' Uniformization
#'
#' @param Q rate matrix.
#' @param t_pow the t in exp(tQ)
#' @param v initial distribution, must be sparse of class dgCMatrix
#' @param eps desired error tolerance
#' @param K precision level. If provided, eps is ignored.
#' @returns v(t) = v %*% exp(tQ) up to tolerance eps.
#' @export
unif_vtexpm = function(Q, K, t_pow=1, v, eps, sparse, verbose=FALSE){
  r = get_max_rate(Q)
  if(missing(K)) K = unif_auto_tune(Q=Q,t_pow=t_pow,eps=eps,r=r)$K
  if(missing(sparse)) sparse = unif_vtexpm_use_sparse(Q)
  if(verbose) cat(sprintf("\t\tunif_vtexpm: N=%d, K=%d, sparse=%s\n",
                          nrow(Q),K,ifelse(sparse,"T","F")))
  nnz = if(sparse) as.double(Matrix::nnzero(Q)) else nrow(Q)^2
  glovars$FLOPS_COUNTER = glovars$FLOPS_COUNTER + K*nrow(v)*nnz
  if(sparse)
    return(unif_vtexpm_sparse(Q, K, t_pow, v, r))
  else
    return(unif_vtexpm_dense(Q, K, t_pow, v, r))
}

###############################################################################
# helpers and wrappers to compiled code
###############################################################################

unif_vtexpm_use_sparse = function(Q){
  (is(Q,"sparseMatrix") && Matrix::nnzero(Q) <= 0.15*(nrow(Q))^2)
}

# vtexpm sparse case
unif_vtexpm_sparse = function(Q, K, t_pow, v, r){
  .Call("R_unif_vtexpm_sparse", PACKAGE = 'pske',
        as(Q,"dgCMatrix"), as.double(t_pow), as.double(r), as(v,"dgCMatrix"),
        nrow(Q), as.integer(K))
}

# vtexpm dense case
unif_vtexpm_dense = function(Q, K, t_pow, v, r){
  .Call("R_unif_vtexpm_dense", PACKAGE = 'pske',
        as(Q,"matrix"), as.double(t_pow), as.double(r), as(v,"matrix"),
        nrow(Q), as.integer(K))
}


