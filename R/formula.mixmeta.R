###
### R routines for the R package mixmeta (c)
#
formula.mixmeta <-
function(x, ...) {
#
################################################################################
# EXTRACTS FORMULA FROM MIXMETA OBJECT (SIMILAR TO lm/glm)
#
  form <- x$formula
  if (!is.null(form)) {
    form <- formula(x$terms)
    environment(form) <- environment(x$formula)
    form
  }
  else formula(x$terms)
}
