###
### R routines for the R package mixmeta (c)
#
model.frame.mixmeta <-
function(formula, ...) {
#
################################################################################
#
  if(is.null(formula$model)) {
    fcall <- formula$call
    fcall$method <- "model.frame"
    env <- environment(formula$terms)
    if(is.null(env)) env <- parent.frame()
    eval(fcall, env)
  }
  else formula$model
}
