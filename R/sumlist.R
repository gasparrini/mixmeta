###
### R routines for the R package mixmeta (c)
#
sumlist <-
function(list) {
#
################################################################################
#
  res <- 0
  for(i in seq(list)) res <- res + list[[i]]
#
  res
}
