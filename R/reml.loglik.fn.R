###
### R routines for the R package mixmeta (c)
#
reml.loglik.fn <-
function(par, Xlist, Zlist, ylist, Slist, nalist, rep, k, q, nall, const,
  bscov, fix) {
#
################################################################################
#
  # COMPUTE Psi FROM PARAMETERS DEPENDING ON STRUCTURE AND PARAMETERIZATION
  Psi <- par2Psi(par,bscov,k,q,fix)
#
  # FIT BY GLS
  Sigmalist <- getSigmalist(Zlist,nalist,Psi,Slist)
  gls <- glsfit(Xlist,ylist,Sigmalist,onlycoef=FALSE)
#
  # RESTRICTED LIKELIHOOD FUNCTION
  # RESIDUAL COMPONENT
  res <- -0.5*(crossprod(gls$invtUy-gls$invtUX%*%gls$coef))
  # DETERMINANT COMPONENTS
  det1 <- -sum(sapply(gls$Ulist,function(U) sum(log(diag(U)))))
  tXWXtot <- sumList(lapply(gls$invtUXlist,crossprod))
  det2 <- -sum(log(diag(chol(tXWXtot))))
#
  # RETURN
  as.numeric(const + det1 + det2 + res)
}
