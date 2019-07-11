###
### R routines for the R package mixmeta (c)
#
ml.newton <-
function(Psi, Xlist, Zlist, ylist, Slist, nalist, rep, k, q, nall, const, bscov,
  fix, control) {
#
################################################################################
#
  # SET INITIAL PARAMETERS AND LIKELIHOOD FUNCTIONS
  # NB: DERIVATIVE ONLY FOR STANDARD MODELS AND UNSTRUCTURED FORM
  par <- unlist(Psi2par(Psi,bscov,k,q,fix))
  fn <- ml.loglik.fn
  gr <- if(is.null(Zlist) && bscov=="unstr") ml.loglik.gr else NULL
#
  # MAXIMIZE
  if(control$showiter) cat("Newton iterations:\n")
  opt <- optim(par=par,fn=fn,gr=gr,Xlist=Xlist,Zlist=Zlist,ylist=ylist,
    Slist=Slist,nalist=nalist,rep=rep,k=k,q=q,nall=nall,const=const,bscov=bscov,
    fix=fix,method="BFGS",control=control$optim,hessian=control$hessian)
#
  # Psi: ESTIMATED BETWEEN-STUDY (CO)VARIANCE MATRIX
  Psi <- par2Psi(opt$par,bscov,k,q,fix)
#
  # RETURN
  list(Psi=Psi,par=opt$par,logLik=opt$value,converged=opt$convergence==0,
    niter=opt$counts[[2]],hessian=opt$hessian)
}
