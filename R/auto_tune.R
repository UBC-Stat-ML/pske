###############################################################################
# skeletoid
###############################################################################

get_max_rate=function(Q){
  if(is(Q,"sparseMatrix"))
    (-min(Q@x)) # much faster than getting the diag first (which is expensive to reconstruct)
  else
    (-min(Matrix::diag(Q)))
}

# get approximation order K that should satisfy error < eps
# uses err <= C*2^{-k}, where C=U(k=0)~(q*t_pow)^2 and U is the upper bound
# in proposition xx, and q = min_i q_ii
skeletoid_auto_tune = function(Q, t_pow, eps=1E-7, max_K=glovars$SKEL_MAX_K, r){
  if(missing(r)) r = get_max_rate(Q)
  eps = max(eps,.Machine$double.eps) # no point in trying to go further than this
  K = as.integer(max(0L,min(ceiling(2*log2(r*t_pow) - log2(eps)), max_K)))
  list(K=K,r=r)
}

# # test
# Q = FiniteBirthDeath$new(S=10,mr=1)$Q
# skeletoid_auto_tune(Q,t_pow=100)
# min(diag(Q))
# skeletoid_auto_tune(matrix(-1),1,eps=1E-15)
# ubound = function(k,barqt = -1){
#   two_pk = 2^(-k)
#   prod_bq_2k = barqt*two_pk
#   exp_bar_two_pk = exp(-prod_bq_2k)
#   barqt*barqt*two_pk*(1 + 0.5*exp_bar_two_pk + 0.25*prod_bq_2k*prod_bq_2k*exp_bar_two_pk)
# }
# curve({ubound(x)},from = 0,to=20)

###############################################################################
# uniformization
###############################################################################

unif_auto_tune = function(Q, t_pow, eps=1E-7, max_K=glovars$UNIF_MAX_K, r){
  if(missing(r)) r = get_max_rate(Q)
  if(eps>1) return(list(K=0L,r=r))
  eps = max(eps,10*.Machine$double.eps) # no point in trying to go further than this (warning: don't set it to .Machine$double.eps otherwise we get Inf in qpois)
  lambda=r*t_pow
  K=as.integer(min(max_K,qpois(p = eps,lambda = lambda,lower.tail = FALSE)))
  return(list(K=K,r=r))
}

# # test
# n=100L; t_pow=20; eps=1-1E-1; pnz=6/n; max_K=60L
# Q=FiniteRndSparse$new(S=n,pnz_per_row=pnz)$Q # sample a new matrix
# K=unif_auto_tune(Q,t_pow,eps,max_K = max_K)
# K
# stopifnot(K==max_K || ppois(K,lambda = -min(Matrix::diag(Q))*t_pow,lower.tail = FALSE) <= eps)
