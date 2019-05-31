\name{coef.mixmeta}
\alias{coef.mixmeta}
\alias{vcov.mixmeta}


\title{ Extract Coefficients and (Co)Variance Matrix from mixmeta Objects }

\description{
These method functions return the estimated fixed-effects coefficients and their (co)variance matrix for fitted meta-analytical models represented in objects of class \code{"mixmeta"}.
}

\usage{
\method{coef}{mixmeta}(object, format=c("vector","matrix"), \dots)

\method{vcov}{mixmeta}(object, \dots)
}

\arguments{
  \item{object }{ an object of class \code{"mixmeta"}.}
  \item{format }{ format of the returned object.}
  \item{\dots }{ further arguments passed to or from other methods.}
}

\value{
For \code{coef}, by default a vector (default)with the estimated fixed-effects coefficients. For multivariate models, a matrix can also be returned.

For \code{vcov}, the (co)variance matrix of the estimated fixed-effects coefficients.
}

\author{Antonio Gasparrini <\email{antonio.gasparrini@lshtm.ac.uk}>}

\seealso{
See \code{\link{mixmeta-package}} for an overview of the package and modelling framework.
}

\examples{
# RUN THE MODEL
model <- mixmeta(cbind(PD,AL) ~ pubyear, S=berkey98[5:7], data=berkey98)

# COEFFICIENTS
model$coef
coef(model)
coef(model, format="matrix")
summary(model)$coef

# (CO)VARIANCE MATRIX
vcov(model)
}

\keyword{models}
\keyword{regression}
\keyword{multivariate}
\keyword{methods}
