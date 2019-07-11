###
### R routines for the R package mixmeta (c)
#
getSigmalist <-
function(Zlist, nalist, Psi, Slist)  {
#
################################################################################
# FUNCTION TO DEFINE THE LIST OF MARGINAL (CO)VARIANCE MATRICES
#
  # IF NO Psi, SIMPLY RETURN Slist
  if(is.null(Psi)) return(Slist)
  Psi <- getList(Psi)
#
  # IF NO Zlist, USE INFO IN nalist TO EXCLUDE MISSING
  if(is.null(Zlist)) return(mapply(function(S,na)
    S+Psi[[1L]][!na,!na,drop=FALSE],Slist,nalist,SIMPLIFY=FALSE))
#
  # OTHERWISE, IN EACH LEVEL MULTIPLY BY Z, BLOCK-DIAG, ADD BY LEVEL AND S
  Psi <- getList(Psi)
  lapply(seq_along(Zlist), function(i)
    sumList(lapply(seq_along(Psi),function(j)
      bdiagMat(lapply(Zlist[[i]][[j]],function(x)
        x%*%Psi[[j]]%*%t(x)))))+Slist[[i]])
}
