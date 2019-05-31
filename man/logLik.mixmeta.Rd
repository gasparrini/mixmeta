\name{logLik.mixmeta}
\alias{logLik.mixmeta}


\title{ Extract Log-Likelihood from mixmeta Objects }

\description{
This method function returns the (restricted) log-likelihood for fitted meta-analytical models represented in objects of class \code{"mixmeta"}.
}

\usage{
\method{logLik}{mixmeta}(object, \dots)
}

\arguments{
  \item{object }{ an object of class \code{"mixmeta"}.}
  \item{\dots }{ further arguments passed to or from other methods.}
}

\value{
A numeric scalar of class \code{"logLik"} with attributes, providing the (restricted) log likelihood of the model. Attributes correspond to the component \code{df} of \code{mixmeta} objects, namely the following scalars: \code{nall} (number of observations used for estimation, excluding missing values), \code{nobs} (equal to \code{nall}, minus the number of fixed-effects coefficients for REML models, \code{fixed} (number of estimated fixed-effects coefficients), \code{random} (number of estimated (co)variance terms).
}

\references{
Gasparrini A, Armstrong B, Kenward MG (2012). Multivariate meta-analysis for non-linear and other multi-parameter associations. \emph{Statistics in Medicine}. \bold{31}(29):3821--3839. [Freely available \href{http://www.ag-myresearch.com/2012_gasparrini_statmed.html}{\bold{here}}].
}

\author{Antonio Gasparrini <\email{antonio.gasparrini@lshtm.ac.uk}>}

\note{
This functions is called by \code{\link{AIC}} and \code{\link{BIC}} for computing the Akaike and Bayesian information criteria.
}

\seealso{
See the default method \code{\link{logLik}}. See \code{\link{mixmeta-package}} for an overview of the package and modelling framework.
}

\examples{
# RUN THE MODEL
model <- mixmeta(cbind(PD,AL)~pubyear,S=berkey98[5:7],data=berkey98)

# LOG-LIKELIHOOD
ll <- logLik(model)
ll
attributes(ll)

# AIC and BIC
AIC(model)
BIC(model)
}

\keyword{models}
\keyword{regression}
\keyword{multivariate}
\keyword{methods}
