###############################################################################
# skeletoid: solution to the evolution equation (or chemical master equation)
###############################################################################

#' Skeletoid
#'
#' @param Q rate matrix.
#' @param t_pow the t in exp(tQ)
#' @param v initial distribution, must be sparse of class dgCMatrix
#' @param eps desired error tolerance
#' @param K precision level. If provided, eps is ignored.
#' @returns v(t) = v %*% exp(tQ) up to tolerance eps.
#' @export
skeletoid_vtexpm = function(Q, K, t_pow, v, eps, verbose = FALSE,sparse,
                            n3=glovars$N3, part_K){
  N = nrow(Q); nr=nrow(v)
  if(missing(K)) K = skeletoid_auto_tune(Q = Q,t_pow = t_pow,eps = eps)$K
  if(missing(part_K)) part_K = get_opt_partition(K=K,N=N,n3=n3,nr=nr)
  sp_cost = nr*Matrix::nnzero(Q)*(2^K)
  if(missing(sparse)) sparse=skeletoid_vtexpm_use_sparse(Q,K,part_K,sp_cost)
  if(sparse){ # use only vector*matrix multiplications to preserve sparsity
    k1=0L; k2=K; FLOPS=sp_cost
  }else{
    k1=part_K$k1; k2=part_K$k2; FLOPS=N*N*comp_cost(x=k2,K=K,N=N,n3=1,nr=nr)}
  glovars$FLOPS_COUNTER = glovars$FLOPS_COUNTER + FLOPS
  delta = t_pow*(2^(-K)) # size of increment
  build_S_delta_sp = (sparse || skeletoid_expm_use_sparse(Q)) # should we at least build S_delta using sparse ops?
  if(verbose) cat(sprintf("\t\tskeletoid_vtexpm: N=%d, K=%d, k1=%d, k2=%d, sparse=%s, S_delta sp.=%s\n",
                          N,K,k1,k2,ifelse(sparse,"T","F"),ifelse(build_S_delta_sp,"T","F")))
  if(build_S_delta_sp)
    return(skeletoid_vtexpm_sparse(Q, delta, v, K, k1, k2,full_sp = sparse))
  else
    return(skeletoid_vtexpm_dense(as(Q,"matrix"), delta, as(v,"matrix"), k1, k2))
}

###############################################################################
# skeletoid: helpers and wrappers to compiled code
###############################################################################

# should we use sparse operations?
SK_VTEXPM_SWITCH_SP_DENSE_THRSH=0.15 # from: test_switch_sparse_dense.R
skeletoid_vtexpm_use_sparse = function(Q,K,part_K,sp_cost,
                                       cr_thresh=SK_VTEXPM_SWITCH_SP_DENSE_THRSH){
  if(!is(Q,"sparseMatrix")) return(FALSE)
  if(missing(sp_cost)) sp_cost = Matrix::nnzero(Q)*(2^K)
  dense_cost = part_K$cost
  return(sp_cost/dense_cost <= cr_thresh)
}

# dense case
skeletoid_vtexpm_dense = function(Q, delta, v, k1, k2){
  .Call("R_skeletoid_vtexpm_dense", PACKAGE = 'pske',
        as.double(Q), as.double(delta), as.double(v),
        nrow(v), ncol(Q), as.integer(k1), as.integer(k2))
}

# sparse case
skeletoid_vtexpm_sparse = function(Q, delta, v, K, k1, k2, full_sp){
  .Call('R_skeletoid_vtexpm_sparse', PACKAGE = 'pske',
        as(Q,"dgCMatrix"), as.double(delta), as(v,"dgCMatrix"), as.integer(K),
        as.integer(k1), as.integer(k2), as.logical(full_sp))
}

#######################################
# optimally distributing work between matrix-matrix and vector-matrix mult
# n3: "effective" number of O(N^3) operations inside squaring loop (the one for k1)
#     Accounts for the fact that multiplying 2 (n \times n) matrices together
#     is faster than doing n row-vector * matrix ops, 1 for each row of the 1st mat.
# NOTE: see test_skeletoid_vtexpm_pick_n3 on choosing n3
#######################################

# computing cost up to scale, actual cost is N^2 * [n3*(K-x)*N + nr*(2^x)]
comp_cost = function(x,K,N, n3,nr){ n3*(K-x)*N + nr*(2^x) }

# optimally partition K to minimize computing cost
# derivation:
# 0 = d/dx Cost(x) = -n3*N + nr*log(2)2^x => (n3*N)/(nr*log(2))=2^x 
# => x* = log2(n3*N/nr) - log2*(log(2))
NEG_LOG2_LOG2 = 0.5287663729448977 # = -log2(log(2))
get_opt_partition = function(K,N, n3=glovars$N3,nr=1){
  x = log2(n3*N/nr) + NEG_LOG2_LOG2 # real valued candidate
  if(x<=0){
    k2=0L
    n_ops = N*N*comp_cost(0,K,N,n3,nr)
  } else if(x>K){
    k2 = K
    n_ops = N*N*comp_cost(K,K,N,n3,nr)
  } else{
    trunc_x=c(ceiling(x),floor(x)) # truncate both ways
    cost_x = comp_cost(x=trunc_x, K=K,N=N,n3=n3,nr=nr) # get both costs
    ind = which.min(cost_x) # find which is better
    k2 = as.integer(trunc_x[ind]) # force integer
    n_ops=N*N*cost_x[ind]
  }
  return(list(k1 = K-k2, k2 = k2, cost = n_ops))
}

# # test
# get_opt_partition(10,50,0.15)

#######################################
# find smallest K such that cost <= given cost
# DENSE CASE:
# for given K, real-valued optimum occurs at x* = log2(n3*N/(nr*log(2))).
# Optimal cost is then:
# = N^2*[ n3*(K-log2(n3*N/(nr*log(2))))*N + 2^{log2(n3*N/(nr*log(2)))} ]
# = N^2*[ n3*(K-log2(n3*N/(nr*log(2))))*N + n3*N/(nr*log(2)) ]
# = n3*(N^3)*(K-log2(n3*N/(nr*log(2))) + 1/(nr*log(2)))
# Hence, (C/n3)N^{-3} = K_match-log2(n3*N/(nr*log(2))) + 1/(nr*log(2))
# => K_match = (C/n3)N^{-3} + x* - 1/(nr*log(2))
#
# SPARSE CASE: let c_r = SK_VTEXPM_SWITCH_SP_DENSE_THRSH (defined above)
# cost_u*c_r = cost_sk_sparse <=> nr*nnz*K_u*c_r=nr*nnz*2^K_s <=> K_s=as.integer(log2(K_u*c_r))
#######################################

INV_LOG_2 = 1.442695040888963 # = 1/log(2)
get_K_match_cost = function(w,N,sparse,nnz,K_u,un_use_sparse,nr,n3=glovars$N3,
                            c_r=SK_VTEXPM_SWITCH_SP_DENSE_THRSH){
  if(w<=2) return(list(K_s=0L,sparse=FALSE))
  # compute the 3 dense options
  vec_mat_tr_cost = (w/N)/N
  mat_mat_tr_cost = vec_mat_tr_cost/(n3*N)
  x_opt = log2(n3*N/nr) + NEG_LOG2_LOG2
  # note: as.integer does implicit "floor". We need it because otherwise we violate cost restriction
  K_vec=as.integer(c(
    mat_mat_tr_cost + x_opt - INV_LOG_2/nr, # mixed
    mat_mat_tr_cost, # all mat-mat products
    log2(vec_mat_tr_cost) # all vector matrix products (dense)
  ))
  if(K_vec[1L]<floor(x_opt)) K_vec[1L]=0L # splitting at x_opt is infeasible
  # sparse option: if given cost w is in sparse ops, convert 1:1. Otherwise, penalize.
  K_sp = if(un_use_sparse) log2(K_u) else log2(K_u*c_r) 
  K_vec = c(as.integer(K_sp),K_vec)
  ind_best = which.max(K_vec)
  return(list(K_s=K_vec[ind_best],sparse=(ind_best==1L)))
}

# # test
# N=500L
# stopifnot(all(
#   get_K_match_cost(w = get_opt_partition(1L,N)$cost,N)==1L,
#   get_K_match_cost(w = get_opt_partition(4L,N)$cost,N)==4L,
#   get_K_match_cost(w = get_opt_partition(18L,N)$cost,N)==18L
# ))
