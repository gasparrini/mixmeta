###
### R routines for the R package mixmeta (c)
#
dropList <-
function(object) {
#
################################################################################
# DROP THE LIST STRUCTURE IF THE LIST HAS ONLY ONE COMPONENT
#
  if(is.list(object) && length(object)==1L) object[[1L]] else object
}
