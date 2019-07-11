\name{summary.mixmeta}
\alias{summary.mixmeta}
\alias{print.mixmeta}
\alias{print.summary.mixmeta}

\title{ Summarizing mixmeta Models }

\description{
Print and summary method functions for fitted meta-analytical models represented in objects of class \code{"mixmeta"}.
}

\usage{
\method{summary}{mixmeta}(object, ci.level=0.95, \dots)

\method{print}{summary.mixmeta}(x, digits=4, report=c("sd","var"), \dots)

\method{print}{mixmeta}(x, digits=4, \dots)
}

\arguments{
  \item{object }{ an object of class \code{"mixmeta"} produced by a call to \code{\link{mixmeta}}.}
  \item{x }{ an object of class \code{"mixmeta"} or \code{"summary.mixmeta"}, produced by calls to \code{\link{mixmeta}} or \code{summary.mixmeta}, respectively.}
  \item{ci.level }{ a numerical value between 0 and 1, specifying the confidence level for the computation of confidence intervals.}
  \item{digits }{ an integer specifying the number of digits to which printed results must be rounded.}
  \item{report }{ if standard deviations (\code{sd}) or variances (\code{var}) must be reported for summarizing the random-effects (co)variance structure.}
  \item{\dots }{ further arguments passed to or from other methods.}
}

\details{
The \code{print} method function for class \code{"mixmeta"} only returns basic information on the fitted model, namely the call, estimated fixed-effects coefficients, dimensions and fit statistics (log-likelihood, AIC, BIC).

The \code{summary} method function computes additional statistics and tests, and produces a list object of class \code{"summary.mixmeta"}. The \code{print} method function for this class shows additional information, such as tables reporting the estimates for the fixed and random-effects parts of the model, Cochran Q test for heterogeneity and I-squared.
}

\value{
The \code{summary} method function for \code{mixmeta} objects produces a list of class \code{"summary.mixmeta"}. The components of the lists are some of those stored in the related \code{mixmeta} object, plus the following:

  \item{coefficients }{ a matrix reporting point estimates, standard errors, z statistics and related p-values of the test, and confidence intervals for the \eqn{kp} fixed-effects coefficients. Note this is different than the component with the same name stored in \code{mixmeta} objects, simply reporting the point estimates (see \code{\link{mixmetaObject}}).}
  \item{AIC }{ the value of the Akaike information criterion for the fitted \code{mixmeta} model, obtained through a call to \code{\link{AIC}}.}
  \item{BIC }{  the value of the Bayesian information criterion for the fitted \code{mixmeta} model, obtained through a call to \code{\link{BIC}}.}
  \item{corFixed }{ the \eqn{kp \times kp}{kp x kp} correlation matrix of the fixed-effects coefficients, obtained from the (co)variance matrix \code{vcov} (see \code{\link{mixmetaObject}} and \code{\link[=vcov.mixmeta]{vcov}}).}
  \item{corRandom }{ the \eqn{kq \times kq}{kq x kq} correlation matrix of the random effects, obtained from the random-effects (co)variance matrix \code{Psi}, or a list of multiple matrices for multilevel models. See \code{\link{mixmetaObject}}.}
  \item{qstat }{ results from the Cochran Q test for heterogeneity, namely a list corresponding to a \code{qtest.mixmeta} object without its class, obtained through \code{\link[=qtest.mixmeta]{qtest}}.}
  \item{i2stat }{ I-squared statistic for the meta-analytical model.}
  \item{ci.level }{ the confidence level used for defining the confidence intervals for the estimates of the fixed-effects coefficients.}

As usual, the \code{print} method functions for classes \code{"mixmeta"} and \code{"summary.mixmeta"} do not return any value.
}

\references{
Sera F, Gasparrini A (2019). An extended mixed-effects framework for meta-analysis.\emph{Statistics in Medicine}. In submission.

}

\author{Antonio Gasparrini <\email{antonio.gasparrini@lshtm.ac.uk}> and Francesco Sera <\email{francesco.sera@lshtm.ac.uk}>}

\seealso{
See \code{\link{mixmeta}} and \code{\link{mixmetaObject}}.
}


\examples{
# RUN THE MODEL
model <- mixmeta(cbind(PD,AL) ~ pubyear, S=berkey98[5:7], data=berkey98)

# SIMPLE PRINT
model

# DEFINE DIGITS
print(model, digit=2)

# SUMMARY WITH 80TH CONFIDENCE INTERVALS
summary(model, ci.level=0.80)

# REPORT RANDOM EFFECTS IN TERMS OF VARIANCES (USE print)
print(summary(model), report="var")
}

\keyword{models}
\keyword{regression}
\keyword{multivariate}
\keyword{methods}
