###############################################################################
# elementwise-max of skeletoid and uniformization
# uses the lowest effort possible as "quoted" by autotuners for given tolerance
###############################################################################

###############################################################################
# expm
###############################################################################

#' Point-wise maximum of Skeletoid and Uniformization computed in parallel
#'
#' @param Q rate matrix.
#' @param t_pow the t in exp(tQ)
#' @param eps desired error tolerance
#' @returns exp(tQ) up to tolerance eps in op-1 norm.
#' @export
max_sk_unif_expm = function(Q, t_pow, eps){
  # request an "effort" quote from both methods, select the lowest
  res = skeletoid_auto_tune(Q,t_pow,eps)
  K_s = res$K; r = res$r # get rate of the Poisson process "for free"
  K_u = unif_auto_tune(Q=Q,t_pow=t_pow,eps=eps,r=r)$K
  K_s = K_u = min(K_u,K_s)
  
  # check if solvers would use sparse operations
  sk_use_sparse = skeletoid_expm_use_sparse(Q)
  un_use_sparse = unif_expm_use_sparse(Q)
  Q_list = set_sp_dense(Q,c(sk_use_sparse, un_use_sparse))
  return(max_sk_unif_expm_(Q_list,t_pow,r,nrow(Q),K_u,K_s,un_use_sparse,
                           sk_use_sparse))
}

###############################################################################
# vtexpm
###############################################################################

#' Point-wise maximum of Skeletoid and Uniformization computed in parallel
#'
#' @param Q rate matrix.
#' @param t_pow the t in exp(tQ)
#' @param v initial distribution, must be sparse of class dgCMatrix
#' @param eps desired error tolerance
#' @returns v(t) = v %*% exp(tQ) up to tolerance eps.
#' @export
max_sk_unif_vtexpm = function(Q, t_pow, v, eps, verbose=FALSE){
  sparse = is(Q,"sparseMatrix"); N = nrow(Q)
  nnz = as.double(if(sparse) Matrix::nnzero(Q) else N*N) # cast double to avoid int overflow below
  
  # request an "effort" quote from skeletoid and its cost
  res = skeletoid_auto_tune(Q,t_pow,eps); K_s = res$K; r = res$r
  part_K = get_opt_partition(K=K_s,N=N)
  sk_use_sparse = skeletoid_vtexpm_use_sparse(Q,K_s,part_K)
  cost_sk = if(sk_use_sparse) nnz*(2^K_s) else part_K$cost
  
  # request an "effort" quote from uniformization and its cost
  K_u = unif_auto_tune(Q=Q,t_pow=t_pow,eps=eps,r=r)$K
  cost_u = nnz*K_u
  un_use_sparse = unif_vtexpm_use_sparse(Q)
  
  # "buy" the cheapest, and lower the effort of loser to match winner's cost,
  # so that they finish at the same time
  if(cost_u <= cost_sk){ # uniformization is faster
    res = get_K_match_cost(w=cost_u,N,sparse,nnz,K_u,un_use_sparse) # best of dense and sparse opts
    K_s=res$K_s; sk_use_sparse=res$sparse
    part_K = if(sk_use_sparse) list(k1=0L,k2=K_s) else get_opt_partition(K=K_s,N=N)
  } else{ # skeletoid is faster
    K_u = as.integer(cost_sk/nnz) # implicit floor'ing, since we cant surpass given cost
    # no need to recheck sparse use, rule doesn't depend on K
  }
  sk_use_sparse_full = sk_use_sparse # adapt to C code convention
  sk_use_sparse = (sk_use_sparse_full || skeletoid_expm_use_sparse(Q)) # might just build S_delta sparse and then do all dense
  
  # print if needed
  if(verbose){
    cat(sprintf("\t\tunif_vtexpm: N=%d, K=%d, sparse=%s\n",
                N,K_u,ifelse(un_use_sparse,"T","F")))
    cat(sprintf("\t\tskeletoid_vtexpm: N=%d, K=%d, k1=%d, k2=%d, sparse=%s, S_delta sp.=%s\n",
                N,K_s,part_K$k1,part_K$k2,
                ifelse(sk_use_sparse_full,"T","F"),
                ifelse(sk_use_sparse,"T","F")))
  }
  
  # build input lists and go!
  use_sparse_vec = c(sk_use_sparse, un_use_sparse)
  Q_list = set_sp_dense(Q,use_sparse_vec) # get Q placeholders
  v_list = set_sp_dense(v,use_sparse_vec) # get v placeholders
  return(max_sk_unif_vtexpm_(Q_list,t_pow,r,nrow(v),N,K_u,K_s,part_K,v_list,
                             un_use_sparse,sk_use_sparse,sk_use_sparse_full))
}

###############################################################################
# helpers
###############################################################################

# for Q or v, gives correctly labeled sparse and dense versions, using dummies
# whenever an option is not needed
set_sp_dense = function(A,use_sparse_vec){
  if(any(use_sparse_vec))
    A_sp = A # A must be sparse already
  else
    A_sp = Matrix::sparseMatrix(i=integer(),j=integer(),x=double(),dims = c(0L,0L))
  
  if(any(!use_sparse_vec))
    A_de = as.double(as(A, "matrix")) # this costs 0 if it's already dense
  else
    A_de = matrix(double()) # nobody needs a dense matrix, set to empty matrix

  return(list(A_sp=A_sp,A_de=A_de))
}

# expm caller
max_sk_unif_expm_ = function(Q_list,t_pow,r,n,K_u,K_s,un_use_sparse,
                             sk_use_sparse){
  .Call("R_max_sk_unif_expm",PACKAGE = "ProbSKE",
        as.double(Q_list$A_de),Q_list$A_sp,as.double(t_pow),as.double(r),
        n, as.integer(K_u),as.integer(K_s), as.logical(un_use_sparse),
        as.logical(sk_use_sparse))
}

# vtexpm caller
max_sk_unif_vtexpm_ = function(Q_list,t_pow,r,nr,nc,K_u,K_s,part_K,v_list,
                               un_use_sparse,sk_use_sparse,sk_use_sparse_full){
  .Call("R_max_sk_unif_vtexpm",PACKAGE = "ProbSKE",
        Q_list$A_de,Q_list$A_sp,as.double(t_pow),as.double(r),
        nr,nc, as.integer(K_u),as.integer(K_s), part_K$k1, part_K$k2,
        v_list$A_de,v_list$A_sp, as.logical(un_use_sparse),
        as.logical(sk_use_sparse), as.logical(sk_use_sparse_full))
}