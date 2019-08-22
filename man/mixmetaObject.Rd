\name{mixmetaObject}
\alias{mixmetaObject}


\title{ mixmeta Objects }

\description{
An object returned by the \code{\link{mixmeta}} function, inheriting from class \code{"mixmeta"}, and representing a fitted univariate or multivariate meta-analytical model.
}

\value{
Objects of class \code{"mixmeta"} are lists with defined components. Dimensions of such components may refer to \eqn{k} outcome parameters, \eqn{p} fixed-effects and \eqn{q} random-effects predictors, \eqn{m} groups and \eqn{n} units used for fitting the model (the latter can be different from those originally selected due to missing). Depending on the type of meta-analytical model, the following components can bu included in a legitimate \code{mixmeta} object:

  \item{\code{coefficients} }{ a \eqn{kp}-dimensional vector of the fixed-effects coefficients.}
    \item{vcov }{ estimated \eqn{kp \times kp}{kp x kp} (co)variance matrix of the fixed-effects coefficients.}
  \item{Psi }{ the estimated \eqn{kq \times kq}{kq x kq} random-effects (co)variance matrix, or a list of matrices for multilevel models. Only for random-effects models.}
  \item{residuals }{ a \eqn{n}-dimensional vector (for univariate models) or \eqn{n \times k}{n x k} matrix (for multivariate models) of residuals, that is observed minus fitted values.}
  \item{fitted.values }{ a \eqn{n}-dimensional vector (for univariate models) or \eqn{n \times k}{n x k} matrix (for multivariate models) of fitted mean values.}
  \item{df.residual }{ the residual degrees of freedom.}
  \item{rank }{ the numeric rank of the fixed-effects part of the fitted model.}
  \item{logLik }{ the (restricted) log-likelihood of the fitted model. Set to \code{NA} for non-likelihood models.}
  \item{converged, niter }{ for models with iterative estimation methods, logical scalar indicating if the algorithm eventually converged and number or iterations, respectively.}
  \item{par }{ parameters estimated in the optimization process when using likelihood-based estimators. These correspond to trasformations of entries of the random-effects (co)variance matrix, dependent on chosen \code{\link[=mixmetaCovStruct]{(co)variance structure}}. For multilevel models, the vector includes the parameters of multiple matrices. Returned also for full (R)IGLS optimization, even if not directly used. See also the \code{\link{mixmeta.ml}} for details.}
  \item{hessian }{ Hessian matrix of the estimated parameters in \code{par} above, only returned if \code{hessian=TRUE} in \code{\link{mixmeta.control}}. See the related \code{\link[=ml.loglik.fn]{optimizations algorithms}} for details.}
  \item{dim }{ list with the following components: \code{k} (scalar, number of outcome parameters), \code{n} (scalar, number of units included in estimation, which could be lower than the total number in the presence of missing values), \code{m} (scalar, number of outer-level groups), \code{p} (scalar, number of fixed-effects predictors), \code{q} (scalar or vector,number of random-effects predictors).}
  \item{df }{ list with the following scalar components: \code{nall} (number of observations used for estimation, excluding missing values), \code{nobs} (equal to \code{nall}, minus the number of fixed-effects coefficients in REML models), \code{fixed} (number of estimated fixed-effects coefficients), \code{random} (number of estimated random-effects (co)variance terms).}
  \item{lab }{ list with the following label vectors: \code{k} for the outcome parameters, and \code{p} and \code{q} for the fixed and random-effects predictors, respectively (including intercept). The last one can be a list for multilevel models.}
    \item{S }{ a \eqn{n \times k(k+1)/2}{n x k(k+1)/2} matrix, where each row represents the vectorized entries of the lower triangle of the related within-unit (co)variance error matrix, taken by column. See \code{\link{mixmeta}}.}
  \item{call }{ the function call.}
  \item{formula }{ the formula for the fixed-effects part of the model. See \code{\link{mixmetaFormula}}.}
  \item{model }{ the model frame used for fitting. Reported if \code{model=TRUE} in \code{\link{mixmeta}}. See \code{\link[=model.frame.mixmeta]{model.frame}}.}
  \item{terms }{ the \code{\link[=terms.mixmeta]{terms}} object representing the fixed-effects part of the fitted model.}
  \item{contrasts }{ (where relevant) the contrasts used.}
  \item{xlevels }{ (where relevant) a record of the levels of the factors used in fitting.}
  \item{na.action }{ (where relevant) information returned by \code{\link[=model.frame.mixmeta]{model.frame}} on the special handling of NAs. See info on \code{\link[=na.omit.data.frame.mixmeta]{missing values}}.}
  \item{method }{ the estimation method.}
  \item{random }{ the formula (or list of formulae for multilevel models) for the random-effects part of the model. See \code{\link{mixmetaFormula}}.}
  \item{bscov }{ a string defining the random-effects (co)variance structure in likelihood based models. See \code{\link{model.frame.mixmeta}}.}
  \item{control }{ a list with the values of the control arguments used, as returned by \code{\link{mixmeta.control}}.}
}

\section{Methods}{
A number of methods functions are available for \code{mixmeta} objects, most of them common to other regression functions.

Specifically-written method functions are defined for \code{\link[=predict.mixmeta]{predict}} (standard predictions) and \code{\link[=blup.mixmeta]{blup}} (best linear unbiased predictions). The method function \code{\link[=simulate.mixmeta]{simulate}} produces simulated outcomes from a fitted model, while \code{\link[=qtest.mixmeta]{qtest}} performs the Cochran Q test for heterogeneity. Other methods have been produced for \code{\link[=summary.mixmeta]{summary}}, \code{\link[=logLik.mixmeta]{logLik}}, \code{\link[=coef.mixmeta]{coef}}, and \code{\link[=vcov.mixmeta]{vcov}}.

Specific methods are also available for \code{\link[=model.frame.mixmeta]{model.frame}} and \code{\link[=model.matrix.mixmeta]{model.matrix}}. In particular, the former produces the model frame (a data frame with special attributes storing the variables used for fitting) with the additional class \code{"data.frame.mixmeta"}. A method \code{\link[=terms.mixmeta]{terms}} is also available for extracting the terms object (only for fixed-effects or full). Methods \code{\link[=na.omit.data.frame.mixmeta]{na.omit}} and \code{\link[=na.exclude.data.frame.mixmeta]{na.exclude}} for this class are useful for the handling of missing values in \code{mixmeta} objects.

Printing functions for the objects of classes defined above are also provided. \code{anova} methods for performing tests in \code{mixmeta} objects are in development.

All the methods above are visible (exported from the namespace) and documented. In additions, several default method functions for regression are also  applicable to objects of class \code{"mixmeta"}, such as \code{\link{fitted}}, \code{\link{residuals}}, \code{\link{AIC}} and \code{\link{BIC}}, \code{\link{drop1}} and \code{\link{add1}}, or \code{\link{update}}, among others.
}

\references{
Sera F, Gasparrini A (2019). An extended mixed-effects framework for meta-analysis.\emph{Statistics in Medicine}. DOI: 10.1002/sim.8362 [Freely available \href{http://www.ag-myresearch.com/2019_sera_statmed.html}{\bold{here}}].
}

\author{Antonio Gasparrini <\email{antonio.gasparrini@lshtm.ac.uk}> and Francesco Sera <\email{francesco.sera@lshtm.ac.uk}>}

\seealso{
See \code{\link{mixmeta}}. See \code{\link{lm}} or \code{\link{glm}} for standard regression functions. See \code{\link{mixmeta-package}} for an overview of this modelling framework.
}

\examples{
# RUN THE MODEL
model <- mixmeta(cbind(PD,AL)~pubyear, S=berkey98[5:7], data=berkey98)

# INSPECT THE OBJECT
names(model)

# LABELS
model$lab

# FORMULA
model$formula

# CONVERGED?
model$converged
}

\keyword{models}
\keyword{regression}
\keyword{multivariate}

