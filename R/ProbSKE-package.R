#' @useDynLib ProbSKE
#' @import RcppEigen
#' @importFrom Rcpp evalCpp
NULL

# global variables
glovars = new.env()
glovars$FLOPS_COUNTER=0 # FLOPs counter for vtexpm methods
glovars$N3 = 0.1 # (see get_opt_partition)
glovars$UNIF_MAX_K = .Machine$integer.max # max possible value for unif accuracy param
glovars$SKEL_MAX_K = .Machine$integer.max # max possible value for skel accuracy param