###
### R routines for the R package mixmeta (c)
#
logLik.mixmeta <-
function(object, ...) {
#
################################################################################
#
  val <- object$logLik
  attributes(val) <- object$df
#
  class(val) <- "logLik"
#
  val
}

