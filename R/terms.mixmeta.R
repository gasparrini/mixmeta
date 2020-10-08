###
### R routines for the R package mixmeta (c)
#
terms.mixmeta <-
function(x, type="fixed", ...)  {
#
################################################################################
# METHOD TO EXTRACT THE TERMS 
#
  # CHECK TYPE
  type <- match.arg(type, c("fixed","full"))
#
  # SELECT
  tt <- if(type=="full") attr(model.frame(x), "terms") else x$terms
#
  tt
}
