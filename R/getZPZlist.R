###
### R routines for the R package mixmeta (c)
#
getZPZlist <-
function(Zlist, nalist, Psi)  {
#
################################################################################
# FUNCTION TO DEFINE THE LIST OF RANDOM-EFFECTS (CO)VARIANCE MATRICES
#
  # IF NO Psi, SIMPLY RETURN NULL
  if(is.null(Psi)) return(NULL)
#
  # IF NO Zlist, USE INFO IN nalist TO EXCLUDE MISSING
  if(is.null(Zlist)) return(lapply(nalist, function(na) Psi[!na,!na,drop=FALSE]))
#
  # OTHERWISE, IN EACH LEVEL MULTIPLY BY Z, BLOCK-DIAG, ADD BY LEVEL
  Psi <- getList(Psi)
  lapply(seq_along(Zlist), function(i)
    sumList(lapply(seq_along(Psi),function(j)
      bdiagMat(lapply(Zlist[[i]][[j]],function(x)
        x%*%Psi[[j]]%*%t(x))))))
}
