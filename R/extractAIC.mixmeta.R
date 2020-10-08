###
### R routines for the R package mixmeta (c)
#
extractAIC.mixmeta <-
function(fit, scale, k = 2, ...) {
#
################################################################################
# EXTRACTS THE NUMBER OF OBSERVATIONS USED FOR FITTING. USED BY BIC
#
  c(fit$df$df,AIC(fit))
}
