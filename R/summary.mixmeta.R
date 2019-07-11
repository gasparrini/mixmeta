###
### R routines for the R package mixmeta (c)
#
summary.mixmeta <-
function(object, ci.level=0.95, ...) {
#
################################################################################
#
  # CHECK ci.level
  if(ci.level<=0||ci.level>=1) stop("'ci.level' must be within 0 and 1")
#
  # EXTRACT QUANTITIES
  ind <- as.numeric(matrix(seq(object$coefficients),object$dim$p,object$dim$k,
    byrow=TRUE))
  coef <- object$coefficients[ind]
  vcov <- object$vcov[ind,ind,drop=FALSE]
  dim <- object$dim
  Psi <- object$Psi
  lab <- object$lab
#
###########################################################################
# FIXED EFFECTS
#
  # COMPUTE STATISTICS FOR FIXED EFFECTS
  coef.se <- sqrt(diag(vcov))
  zval <- coef/coef.se
  zvalci <- qnorm((1-ci.level)/2,lower.tail=FALSE)
  pvalue <- 2*(1-pnorm(abs(zval)))
  ci.lb <- coef-zvalci*coef.se
  ci.ub <- coef+zvalci*coef.se
  cilab <- paste(signif(ci.level,2)*100,"%ci.",c("lb","ub"),sep="")
#
  # GENERATE TABLE AS MATRIX
  tabfixed <- cbind(coef,coef.se,zval,pvalue,ci.lb,ci.ub)
   dimnames(tabfixed) <- list(if(dim$k>1L) names(coef) else lab$p,
     c("Estimate","Std. Error","z","Pr(>|z|)",cilab))
#
  # CORRELATION MATRIX OF FIXED EFFECTS
  corFixed <- vcov/outer(coef.se,coef.se)
#
###########################################################################
# RANDOM EFFECTS
#
  corRandom <- if(object$method!="fixed") {
    Psi <- getList(Psi)
    cor <- lapply(Psi,function(x) {
      ran.sd <- sqrt(diag(x))
      x/outer(ran.sd,ran.sd)
    })
    if(length(cor)==1L) cor[[1L]] else cor
  } else NULL
#
###########################################################################
# QTEST AND I2 STATISTICS
#
  qstat <- unclass(qtest(object))
  i2stat <- pmax((qstat$Q-qstat$df)/qstat$Q*100,0)
#
###########################################################################
#
  # DEFINE THE LIST
  keep <- match(c("vcov","Psi","bscov","random","df.res","rank","logLik",
    "converged","niter","method","dim","df","lab","na.action","call","terms"),
    names(object),0L)
  out <- c(list(coefficients=tabfixed),object[keep],list(AIC=AIC(object),
    BIC=BIC(object),corFixed=corFixed,corRandom=corRandom,qstat=qstat,
    i2stat=i2stat,ci.level=ci.level))
#
  class(out) <- "summary.mixmeta"
#
  out
}
