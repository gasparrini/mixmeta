###
### R routines for the R package mixmeta (c)
#
ml.loglik.fn <-
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
  # LIKELIHOOD FUNCTION
  # RESIDUAL COMPONENT
  res <- -0.5*(crossprod(gls$invtUy-gls$invtUX%*%gls$coef))
  # DETERMINANT COMPONENT
  det <- -sum(sapply(gls$Ulist,function(U) sum(log(diag(U)))))
#
  # RETURN
  as.numeric(const + det + res)
}
