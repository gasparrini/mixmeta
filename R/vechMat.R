###
### R routines for the R package mixmeta (c)
#
vechMat <-
function(mat, diag=TRUE) {
#
################################################################################
#
  if(!is.matrix(mat)) mat <- as.matrix(mat)
  if(diff(dim(mat))!=0) stop("Non-square matrix")
#
  mat[lower.tri(mat,diag=diag)]
}
