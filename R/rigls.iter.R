###
### R routines for the R package mixmeta (c)
#
rigls.iter <-
function(Psi, Qlist, Xlist, Zlist, ylist, Slist, nalist, rep, k, q, bscov,
  fix, control) {
#
################################################################################
#
  # IF ALL fixed, RETURN THE SAME Psi
  if(all(fix <- bscov%in%"fixed")) return(Psi)
#
  # FIT BY GLS AND DERIVE THE INVERTED Sigmalist
  Sigmalist <- getSigmalist(Zlist,nalist,Psi,Slist)
  gls <- glsfit(Xlist,ylist,Sigmalist,onlycoef=FALSE)
  invSigmalist <- lapply(gls$invUlist,tcrossprod)
  #
  # RESPONSE VECTORS WITH RESIDUALS MINUS THE WITHIN (CO)VARIANCE, PLUS REML TERM
  invtXinvSigmaX <- solve(crossprod(gls$invtUX))
  flist <- mapply(function(y,S,X) tcrossprod(y-X%*%gls$coef)-S+
      X%*%invtXinvSigmaX%*%t(X),ylist,Slist,Xlist,SIMPLIFY=FALSE)
#
  # DEFINE THE COMPONENTS (GOLDSTEIN COMP STAT & DATA ANAL 1992, PAGE 66)
  Alist <- lapply(seq(Qlist), function(i)
    lapply(seq(Qlist[[1L]]), function(k) Qlist[[i]][[k]]%*%invSigmalist[[i]]))
  Blist <- lapply(seq(flist), function(i) flist[[i]]%*%invSigmalist[[i]])
#
  # COMPUTE THE ELEMENTS
  ind1 <- unlist(lapply(seq(Qlist[[1L]]),":",length(Qlist[[1L]])))
  ind2 <- rep(seq(Qlist[[1L]]),rev(seq(Qlist[[1L]])))
  XtVX <- xpndMat(sapply(seq(ind1), function(k)
    sum(sapply(seq(Qlist), function(i)
      sum(diag(Alist[[i]][[ind1[k]]]%*%Alist[[i]][[ind2[k]]]))))))
  XtVy <- sapply(seq(Qlist[[1L]]), function(k)
    sum(sapply(seq(Qlist), function(i)
      sum(diag(Alist[[i]][[k]]%*%Blist[[i]])))))
#
  # COMPUTE PARAMETERS AND FORM Psi (NB: par NOT TRANSFORMED)
  par <- as.numeric(chol2inv(chol(XtVX)) %*% XtVy)
  Psi <- lapply(lapply(seq_along(q),function(j)
    par[seq(c(0,cumsum(q*k*(q*k+1)/2))[j]+1,cumsum(q*k*(q*k+1)/2)[j])]),xpndMat)
#
  # FORCE POSITIVE-DEFINITENESS
  Psi <- checkPD(Psi,set.negeigen=control$set.negeigen,force=TRUE,error=FALSE)
#
  # DROP THE LIST STRUCTURE IF ONLY ONE COMPONENT
  dropList(Psi)
}
