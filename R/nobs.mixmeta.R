###
### R routines for the R package mixmeta (c)
#
nobs.mixmeta <-
function (object, ...) {
#
################################################################################
# EXTRACTS THE NUMBER OF OBSERVATIONS USED FOR FITTING. USED BY BIC
#
  object$df$nobs
}
