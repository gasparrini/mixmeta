###
### R routines for the R package mixmeta (c)
#
checkPD <-
function(x, set.negeigen=sqrt(.Machine$double.eps), force=TRUE, error=FALSE,
  label="x") {
#
################################################################################
#
  # TRANFORM IN A LIST
  islist <- is.list(x)
  if(!islist) x <- list(x)
#
  # CHECK POSITIVE-DEFINITENESS
  x <- lapply(x, function(mat) {
    # COMPUTE EIGENVALUES
    eig <- eigen(mat)
    # CHECK IF POSITIVE, AND RETURN ERROR IF REQUIRED
    if(any(ind <- eig$values<0) && error)
      stop(paste("Problems with positive-definiteness in '",label,"'. ",sep=""))
    # FORCE POSITIVENESS, IF REQUIRED
    if(any(ind) && force) {
      eig$values[ind] <- set.negeigen
      mat <- eig$vectors %*% diag(eig$values,ncol(mat)) %*% t(eig$vectors)
    }
    return(mat)
  })
#
  # TRANSFORM IN MATRIX IF ONLY ONE COMPONENT
  if(!islist&&length(x)==1L) x[[1]] else x
}
