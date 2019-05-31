###
### R routines for the R package mixmeta (c)
#
mixmeta <-
function(formula, S, data, random, method="reml", bscov="unstr", offset, subset,
  contrasts=NULL, na.action, model=TRUE, control=list()) {
#
################################################################################
#  RESET SOME ARGUMENTS AND FORMULAE FOR FIXED AND RANDOM PARTS
#
  # RESET FORMULA (FIXED), PRESERVING THE ENVIRONMENT
  if (missing(data)) data <- parent.frame()
  if(class(eval(substitute(formula),data))!="formula") {
    formula <- as.formula(paste(deparse(substitute(formula),width.cutoff=499L),
      "~ 1"),env=parent.frame())
  } else if(length(formula)!=3) stop("'formula' must be a two-sided formula")
#
  # CHECK ESTIMATION METHOD AND PRE-SET control
  method <- match.arg(method,c("fixed","ml","reml","mm","vc","model.frame"))
  control <- do.call("mixmeta.control",control)
#
  # RESET random AND bscov
  random <- getRandom(random,method)
  bscov <- getBSCov(bscov,random,method)
#
################################################################################
# CREATE THE CALL (SPECIAL FORMULA)
#
  # CREATE THE ORIGINAL CALL AND THE MODIFIED VERSION
  call  <- match.call()
  mcall <- match.call(expand.dots=FALSE)
  mn <- match(c("formula", "data", "subset", "weights", "na.action", "offset"),
    names(mcall), 0L)
  mcall <- mcall[c(1L, mn)]
  mcall$drop.unused.levels <- TRUE
  mcall[[1L]] <- as.name("model.frame")
#
  # CREATE THE FULL FORMULA INCLUDING FIXED AND RANDOM TERMS
  mcall$formula <- getFullformula(formula, random)
#
################################################################################
# DERIVE THE MODEL FRAME (SPECIAL HANDLING OF MISSING VALUES) AND TERMS
#
  # NOW KEEP THE MISSING
  mcall$na.action <- "na.pass"
  # CREATE MODEL FRAME WITH ADDITIONAL CLASS
  mf <- eval(mcall, parent.frame())
  class(mf) <- c("data.frame.mixmeta",class(mf))
  # NOW HANDLE THE MISSING
  if(missing(na.action)) na.action <- getOption("na.action")
  if(length(na.action)) mf <- do.call(na.action,list(mf))
  # RETURN mf IF REQUIRED
  if(method=="model.frame") return(mf)
  # SAVE TERMS (ONLY FIXED PART)
  terms <- terms(formula)
#
################################################################################
# GROUPS AND RE-ORDER
#
  # DEFINE GROUPING FACTORS
  groups <- getGroups(random,mf)
#
  # RE-ORDER
  ord <- do.call(order,lapply(seq(ncol(groups)),function(i) groups[,i]))
  groups <- groups[ord,,drop=FALSE]
  mf <- mf[ord,,drop=FALSE]
#
################################################################################
# DERIVE OBJECTS FOR FITTING
#
  # GET DESIGN MATRIX AND RESPONSE (AS MATRIX) FOR FIXED PART
  y <- as.matrix(model.response(mf,"numeric"))
  X <- model.matrix(formula[c(1L,3L)],mf,contrasts)
  offset <- as.vector(model.offset(mf))
  if(!is.null(offset)) {
    if(length(offset)!=NROW(y))
      stop("number of offsets should equal number of observations")
    y <- y - offset
  }
#
  # LIST OF DESIGN MATRICES FOR RANDOM PART (ONLY IF NEEDED)
  Z <- getZ(random,mf,contrasts)
#
  # PRODUCE S AS A MATRIX OF VECTORIZED ENTRIES (IF NEEDED INPUT COVARIANCES)
  # IF PROVIDED THROUGH control, ARRANGED LATER BY getSlist
  S <- eval(call$S,data,parent.frame())
  S <- getS(S,y,attr(mf,"na.action"), if(missing(subset)) NULL else
    eval(call$subset,data,parent.frame()),ord,control$Scor,control$checkPD)
#
################################################################################
# FIT THE MODEL, CALLING mixmeta.fit
#
  # MODEL FIT
  fit <- mixmeta.fit(X,Z,y,S,groups,method,bscov,control)
#
################################################################################
#  COMPLETE
#
  # OFFSET AND REORDER, ALSO DROP DIMENSIONS
  if(!is.null(offset)) fit$fitted.values <- fit$fitted.values+offset
  fit$fitted.values <- fit$fitted.values[order(ord),]
  fit$residuals <- fit$residuals[order(ord),]
  fit$offset <- offset[order(ord)]
  fit$S <- S[order(ord),]
#
  # ADD OTHER COMPONENTS
  fit$call <- call
  fit$formula <- formula
  fit$model <- if(model) mf <- mf[order(ord),,drop=FALSE] else NULL
  fit$terms <- terms
  fit$contrasts <- attr(X,"contrasts")
  fit$xlevels <- .getXlevels(terms(mcall$formula),mf)
  fit$na.action <- attr(mf,"na.action")
  fit$method <- method
  fit$random <- random
  fit$bscov <- bscov
  fit$control <- control
#
  class(fit) <- "mixmeta"
#
  fit
}
