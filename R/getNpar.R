###
### R routines for the R package mixmeta (c)
#
getNpar <-
function(bscov, k, q)  {
#
################################################################################
# FUNCTION TO DETERMINE THE NUMBER OF PARAMETERS FOR EACH RANDOM LEVEL
#
  npar <- sapply(seq_along(bscov), function(i) {
    d <- k*q[i]
    switch(bscov[i],
      unstr = d*(d+1)/2,
      diag = d,
      id = 1,
      cs = 2,
      hcs = d+1,
      ar1 = 2,
      har1 = d+1,
      prop = 1,
      cor = d,
      fixed = 0
    )
  })
#
  names(npar) <- names(bscov)
#
  npar
}
