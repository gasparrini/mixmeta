###
### R routines for the R package mixmeta (c)
#
getZlist <-
function(Z, nay, groups, m, k, q)  {
#
################################################################################
# FUNCTION TO DEFINE THE LIST OF DESIGN MATRICES FOR THE RANDOM PART
#
  # IF NULL, RETURN SO, OTHERWISE TRANSFORM IN LIST
  if(is.null(Z)) return(NULL)
  if(length(q)==1L) Z <- list(Z)
#
  # OTHERWISE, GENERATE THE LIST ACCOUNTING FOR MULTIPLE LEVELS
  # FOR EACH GROUP, A q-LENGTH LIST OF LIST OF MATRICES
  lapply(seq(m),function(i) {
    Zi <- lapply(Z, function(x) x[groups[,1]%in%i,,drop=FALSE])
    gi <- groups[groups[,1]%in%i,,drop=FALSE]
    nayi <- nay[groups[,1]%in%i,,drop=FALSE]
    Zij <- lapply(seq(length(q)),function(j)
      lapply(unique(gi[,j]), function(rep) {
        ind <- gi[,j]%in%rep
        Zind <- Zi[[j]][ind,,drop=FALSE]%x%diag(k)
        Zind[!c(t(nayi[ind,])),,drop=FALSE]
      }))
  })
}
