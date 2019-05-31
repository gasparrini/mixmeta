###
### R routines for the R package mixmeta (c)
#
sumList <-
function(list) {
#
################################################################################
#
  res <- 0
  for(i in seq(list)) res <- res + list[[i]]
#
  res
}
