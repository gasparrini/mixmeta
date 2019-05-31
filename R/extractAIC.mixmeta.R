###
### R routines for the R package mixmeta (c)
#
extractAIC.mixmeta <-
function (object, scale, k = 2, ...) {
#
################################################################################
# EXTRACTS THE NUMBER OF OBSERVATIONS USED FOR FITTING. USED BY BIC
#
   c(object$df$df,AIC(object))
}
