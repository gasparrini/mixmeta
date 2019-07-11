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
  # EXTRACT FULL LIST OF TERMS (FIXED+RANDOM) FROM MODEL FRAME
  tt <- attr(model.frame(x), "terms")
#
  # IF FULL, RETURN
  if(type=="full") return(tt)
#
  # ELSE, IF FIXED TERMS ONLY, REMOVE
  ind <- which(!attr(tt,"term.labels")%in%attr(terms(x$formula),"term.labels"))
  if(length(ind)) tt <- tt[-ind]
#
  tt
}
