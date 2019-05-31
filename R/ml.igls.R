###
### R routines for the R package mixmeta (c)
#
ml.igls <-
function(Psi, Xlist, Zlist, ylist, Slist, nalist, rep, k, q, nall, const, bscov,
  fix, control) {
#
################################################################################
#
  # ONLY FOR UNSTRUCTURED RANDOM FORMS
  if(any(bscov!="unstr")&&control$loglik.iter%in%c("igls","rigls"))
    stop("iterative methods 'igls'-'rigls' only available for bscov='unstr'")
#
  # DESIGN MATRIX MAPPING THE PARAMETERS TO BE ESTIMATED
  # GOLDSTEIN COMP STAT & DATA ANAL 1992, PAGE 65
  Qlist <- getQlist(Zlist,nalist,rep,k,q)
#
  # PRODUCE INITIAL VALUES
  niter <- 0
  converged <- FALSE
  reltol <- control$reltol
#
  # START OPTIMIZATION
  if(control$showiter) cat("IGLS iterations:\n")
  maxiter <- ifelse(control$loglik.iter%in%c("hybrid","newton"),
    control$igls.inititer,control$maxiter)
  while(!converged && niter<maxiter) {
#
    # IF PRINTED
    if(control$showiter) {
      par <- unlist(Psi2par(Psi,bscov,k,q,fix))
      logLik <- ml.loglik.fn(par,Xlist,Zlist,ylist,Slist,nalist,rep,k,q,nall,
        const,bscov,fix)
      cat("iter ",niter,": value ",-logLik,"\n",sep="")
    }
#
    # ITERATION
    old <- unlist(Psi)
    Psi <- igls.iter(Psi,Qlist,Xlist,Zlist,ylist,Slist,nalist,rep,k,q,bscov,
      fix,control)
    niter <- niter+1
#
    # CHECK CONVERGENCE
    converged <- all(abs(unlist(Psi)-old)<reltol*abs(unlist(Psi)+reltol))
    if(control$showiter&&converged) cat("converged\n")
  }
#
  # LOG-LIKELIHOOD
  par <- unlist(Psi2par(Psi,bscov,k,q,fix))
  logLik <- ml.loglik.fn(par,Xlist,Zlist,ylist,Slist,nalist,rep,k,q,nall,const,
    bscov,fix)
#
  # RETURN
  list(Psi=Psi,par=par,logLik=logLik,converged=converged,niter=niter)
}
