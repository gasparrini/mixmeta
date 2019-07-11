###
### R routines for the R package mixmeta (c)
#
getList <-
function(object) {
#
################################################################################
# TRANFORM THE OBJECT IN A LIST
#
  if(is.list(object)) object else list(object)
}
