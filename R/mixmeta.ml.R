###
### R routines for the R package mixmeta (c)
#
mixmeta.ml <-
function(Xlist, Zlist, ylist, Slist, nalist, rep, k, q, nall, bscov, control, ...) {
#
################################################################################
#
  # DEFINE FIXED PARTS, (VERY) INITIAL VALUES, AND LIKELIHOOD CONSTANT
  fix <- getPsifix(control$Psifix,bscov,k,q,control$checkPD)
  Psi <- getInitPsi(control$initPsi,bscov,k,q,fix,control$checkPD)
  const <- -0.5*nall*log(2*pi)
#
  # OPTIMIZE: IGLS AND/OR NEWTON
  lliter <- control$loglik.iter
  if(lliter=="rigls") warning("'igls' used instead than 'rigls")
  if(lliter!="newton") {
    opt <- ml.igls(Psi,Xlist,Zlist,ylist,Slist,nalist,rep,k,q,nall,const,bscov,
      fix,control)
    Psi <- opt$Psi
  }
  if(!lliter%in%c("igls","rigls")) opt <- ml.newton(Psi,Xlist,Zlist,ylist,
    Slist,nalist,rep,k,q,nall,const,bscov,fix,control)
#
  # FIT BY GLS
  Sigmalist <- getSigmalist(Zlist,nalist,opt$Psi,Slist)
  gls <- glsfit(Xlist,ylist,Sigmalist,onlycoef=FALSE)
#
  # COMPUTE (CO)VARIANCE MATRIX OF coef
  qrinvtUX <- qr(gls$invtUX)
  R <- qr.R(qrinvtUX)
  Qty <- qr.qty(qrinvtUX,gls$invtUy)
  vcov <- tcrossprod(backsolve(R,diag(1,ncol(gls$invtUX))))
#
  # COMPUTE RESIDUALS (LATER), FITTED AND RANK
  res <- NULL
  fitted <- lapply(Xlist,"%*%",gls$coef)
  rank <- qrinvtUX$rank
#
  # RETURN
  list(coefficients=gls$coef,vcov=vcov,Psi=opt$Psi,residuals=res,
    fitted.values=fitted,df.residual=nall-rank-length(opt$par),rank=rank,
    logLik=opt$logLik,converged=opt$converged,par=opt$par,hessian=opt$hessian,
    niter=opt$niter)
}
