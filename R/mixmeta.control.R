###
### R routines for the R package mixmeta (c)
#
mixmeta.control <-
function(optim=list(), showiter=FALSE, maxiter=100, initPsi=NULL, Psifix=NULL,
  Scor=NULL, addSlist=NULL, inputna=FALSE, inputvar=10^4, loglik.iter="hybrid",
  igls.inititer=10, hessian=FALSE, vc.adj=TRUE, reltol=sqrt(.Machine$double.eps),
  checkPD=NULL, set.negeigen=sqrt(.Machine$double.eps)) {
#
################################################################################
# SET CONTROL PARAMETERS FOR MODEL FITTING, WITH SPECIFIC DEFAULT VALUES
#
  # OPTIM:
  optim <- modifyList(list(fnscale=-1,maxit=maxiter,reltol=reltol),optim)
  if(showiter) {
    optim$trace <- 6
    optim$REPORT <- 1
  }
#
  # ITERATION METHOD FOR LOG-LIKELIHOOD MODELS
  loglik.iter <- match.arg(loglik.iter,c("hybrid","newton","igls","rigls"))
  if(igls.inititer<=0L) igls.inititer <- 0
#
  # RETURN
	list(optim=optim,showiter=showiter,maxiter=maxiter,hessian=hessian,
    initPsi=initPsi,Psifix=Psifix,Scor=Scor,addSlist=addSlist,inputna=inputna,
	  inputvar=inputvar,loglik.iter=loglik.iter,igls.inititer=igls.inititer,
	  vc.adj=vc.adj,reltol=reltol,checkPD=checkPD,set.negeigen=set.negeigen)
}
