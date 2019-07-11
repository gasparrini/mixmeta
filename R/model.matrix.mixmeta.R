###
### R routines for the R package mixmeta (c)
#
model.matrix.mixmeta <-
function(object, ...) {
#
################################################################################
#
  # DROP THE RANDOM TERMS AND CONTRASTS
  tt <- object$terms
  contr <- getContrXlev(tt, object$contrasts)
#
  # RUN THE DEFAULT METHOD
  model.matrix(tt, data=model.frame(object), contrasts.arg=contr)
}
