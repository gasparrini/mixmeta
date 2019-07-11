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
  x <- getList(x)
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
  # DROP THE LIST STRUCTURE IF ONLY ONE COMPONENT
  dropList(x)
}
