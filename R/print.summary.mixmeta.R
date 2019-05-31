###
### R routines for the R package mixmeta (c)
#
print.summary.mixmeta <-
function(x, digits=4, report=c("sd","var"), ...) {
#
################################################################################
#
  # CREATE USEFUL OBJECTS
  methodname <- c("reml","ml","fixed","mm","vc")
  methodlabel <- c("REML","ML","Fixed","Method of moments",
    "Variance components")
  bscovname <- c("unstr","diag","id","cs","hcs","ar1","har1","prop","cor","fixed")
  bscovlabel <- c("General positive-definite","Diagonal",
    "Multiple of identity","Compound symmetry","Heterogeneous compound symmetry",
    "Autoregressive of first order","Heterogeneous autoregressive of first order",
    "Proportional to fixed matrix","Fixed correlation","Fixed")
  int <- attr(x$terms,"intercept")==1L
#
################################################################################
# HEADING AND SUBHEADING
#
  # HEADING
  cat("Call:  ",paste(deparse(x$call),sep="\n",collapse="\n"),"\n\n",sep="")
#
  # SUB-HEADING
  cat(if(x$dim$k==1L) "Uni" else "Multi","variate ",
    ifelse(is.null(x$random),"","extended "),
    ifelse(x$method=="fixed","fixed","random"),"-effects meta-",
    ifelse(x$dim$p-int>0,"regression","analysis"),"\n",sep="")
  # CHECK LATER FOR CHOICE META-ANALYSIS OR METAREGRESSION
  cat("Dimension: ",x$dim$k,"\n",sep="")
  if(x$method!="fixed") {
    cat("Estimation method: ",methodlabel[which(x$method==methodname)],
      "\n",sep="")
  }
  cat("\n")
#
###########################################################################
# FIXED EFFECTS ESTIMATES
#
  cat("Fixed-effects coefficients","\n",sep="")
#
  # COMPUTE SIGNIFICANCE
  signif <- symnum(x$coefficients[,"Pr(>|z|)"],corr=FALSE,na=FALSE,cutpoints=c(0,
    0.001,0.01,0.05,0.1,1),symbols=c("***","**","*","."," "))
#
  # PRODUCE THE TABLE
  tabletot <- formatC(x$coefficients,digits=digits,format="f")
  tabletot <- cbind(tabletot,signif)
  colnames(tabletot)[7] <- ""
#
  # FOR UNIVARIATE MODELS OR SIMPLE META-ANALYSIS
  if(x$dim$k==1L || x$dim$p==1L) {
    print(tabletot,quote=FALSE,right=TRUE,print.gap=2)
  } else {
    p <- x$dim$p
    for(i in seq(x$dim$k)) {
      ind <- seq((i-1)*p+1,(i-1)*p+p)
      table <- tabletot[ind,,drop=FALSE]
      rownames(table) <- x$lab$p
      cat(" ",x$lab$k[i],":","\n")
      print(table,quote=FALSE,right=TRUE,print.gap=2)
    }
  }
  cat("---\nSignif. codes: ",attr(signif,"legend"),"\n\n")
#
###########################################################################
# RANDOM COMPONENTS
#
  if(!x$method=="fixed") {
    cat("Random-effects (co)variance components", "\n", sep="")
    Psi <- getList(x$Psi)
    random <- getList(x$random)
    report <- match.arg(report, c("sd","var"))
#
    # LOOP
    for(j in seq(x$bscov)) {
      if(!is.null(random[[j]]))
        cat(" Formula: ",deparse(random[[j]]),"\n",sep="")
      cat(" Structure: ",bscovlabel[which(x$bscov[j]==bscovname)],"\n",sep="")
#
      # STANDARD DEVIATIONS
      dd <- if(report=="sd") cbind("Std. Dev"=sqrt(diag(Psi[[j]]))) else
        cbind("Var"=diag(Psi[[j]]))
      dd <- formatC(dd,digits=digits,format="f")
      if(length(dd)==1L) rownames(dd) <- ""
#
      # CORRELATIONS
      if(nrow(Psi[[j]])>1L) {
        corRan <- if(length(x$bscov)>1L) x$corRandom[[j]] else x$corRandom
        corRan[upper.tri(corRan,diag=TRUE)] <- NA
        dimnames(corRan) <- NULL
        corRan <- format(corRan[-1,-ncol(corRan),drop=FALSE],digits=digits,
          format="f")
        corlab <- rownames(Psi[[j]])
        corRan <- rbind(corlab[-length(corlab)],corRan)
        colnames(corRan) <- c("Corr",rep("",length(corlab)-2))
        corRan[grep("NA",corRan)] <- ""
      } else corRan <- NULL
#
      # PRODUCE THE TABLE
      print(cbind(dd,corRan),quote=FALSE,right=TRUE,na.print="",print.gap=2)
      cat("\n")
    }
  }
#
###########################################################################
# OVERALL QTEST AND I-SQUARE
#
  Q <- formatC(x$qstat$Q,digits=digits,format="f")
  pvalue <- formatC(x$qstat$pvalue,digits=digits,format="f")
  i2 <- formatC(x$i2stat,digits=1,format="f")
  cat(if(x$qstat$k==1) "Uni" else "Multi","variate ","Cochran Q-test for ",
    if(x$qstat$residual) "residual ", "heterogeneity:","\n",sep="")
  cat("Q = ",Q[1]," (df = ",x$qstat$df[1],"), p-value = ",pvalue[1],"\n",sep="")
  cat("I-square statistic = ",i2[1],"%","\n\n",sep="")

###########################################################################
# FIT STATS
#
  cat(x$dim$n," units, ",x$dim$k," outcome",ifelse(x$dim$k>1L,"s, ",", "),
    x$df$nall," observations, ",x$df$fixed," fixed and ",
    x$df$random," random-effects parameters","\n",sep="")
  if(na <- length(x$na.action)) cat("(",na," unit",ifelse(na>1L,"s",""),
    " removed due to missingness",")\n",sep="")
  if(!x$method%in%c("mm","vc")) {
    table <- c(x$logLik,x$AIC,x$BIC)
    names(table) <- c("logLik","AIC","BIC")
    table <- formatC(table,digits=digits,format="f")
    print(table,quote=FALSE,right=TRUE,print.gap=2)
  }
  cat("\n")
}
