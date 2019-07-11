###
### R routines for the R package mixmeta (c)
#
drop1.mixmeta <-
function(object, ...) {
#
################################################################################
# COMPUTE TABLE OF CHANGE IN FIT WHEN DROPPING TERMS. ONLY FOR FIXED AND ML
#
  # CHECK ESTIMATION METHOD
  if(!object$method %in% c("fixed","ml")) 
    stop("Fit only comparable in models with method='fixed' or method='ml'")
#
  # RUN GENERIC
  NextMethod(generic="drop1", ...)
}
