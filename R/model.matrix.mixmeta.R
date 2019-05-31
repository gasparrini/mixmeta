###
### R routines for the R package mixmeta (c)
#
model.matrix.mixmeta <-
function(object, ...) {
#
################################################################################
#
  model.matrix(object$formula[c(1,3)],data=model.frame(object),
    contrasts.arg=object$contrasts)
}
