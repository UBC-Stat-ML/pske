##############################################################################
# reset and retrieve the ops counter
##############################################################################

#' @export
reset_ops_counter = function(){
    glovars$FLOPS_COUNTER = 0
    return(NULL)
}

#' @export
get_ops_counter = function(){
    return(glovars$FLOPS_COUNTER)
}

