###
### R routines for the R package mixmeta (c)
#
coef.mixmeta <-
function(object, format=c("vector","matrix"), ...) {
#
################################################################################
#
  coef <- object$coefficients
  format <- match.arg(format,c("vector","matrix"))
  if(format=="matrix" && object$dim$k>1L) coef <- matrix(coef,ncol=object$dim$k,
    dimnames=list(object$lab$p,object$lab$k))
#
  coef
}
