###
### R routines for the R package mixmeta (c)
#
print.mixmeta <-
function(x, digits=4, ...) {
#
################################################################################
# HEADING AND SUB-HEADING
#
  # HEADING
  cat("Call:  ",paste(deparse(x$call),sep="\n",collapse="\n"),"\n\n",sep="")
#
  # SUB-HEADING
  cat("Fixed-effects coefficients:","\n",sep="")
  coef <- if(x$dim$k==1L||x$dim$p==1L) x$coefficients else
    matrix(x$coefficients,x$dim$p,x$dim$k,byrow=TRUE,
      dimnames=list(x$lab$p,x$lab$k))
  table <- formatC(coef,digits=digits,format="f")
  print(table,quote=FALSE,right=TRUE,print.gap=2)
  cat("\n")
#
################################################################################
# FIT STATS
#
  cat(x$dim$n," units, ",x$dim$k," outcome",ifelse(x$dim$k>1L,"s, ",", "),
    x$df$nall," observations, ",x$df$fixed," fixed and ",
    x$df$random," random-effects parameters","\n",sep="")
  if(na <- length(x$na.action)) cat(" (",na," unit",ifelse(na>1L,"s",""),
    " removed due to missingness",")\n",sep="")
  if(!x$method%in%c("mm","vc")) {
    table <- c(x$logLik,AIC(x),BIC(x))
    names(table) <- c("logLik","AIC","BIC")
    table <- formatC(table,digits=digits,format="f")
    print(table,quote=FALSE,right=TRUE,print.gap=2)
  }
  cat("\n")
}
