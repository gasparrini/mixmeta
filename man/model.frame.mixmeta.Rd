\name{model.frame.mixmeta}
\alias{model.frame.mixmeta}
\alias{model.matrix.mixmeta}


\title{ Extract Model Frame and Design Matrix from mixmeta Objects }

\description{
These method functions return the model frame and design matrix for meta-analytical models represented in objects of class \code{"mixmeta"}.
}

\usage{
\method{model.frame}{mixmeta}(formula, \dots)

\method{model.matrix}{mixmeta}(object, \dots)
}

\arguments{
  \item{object, formula }{ an object of class \code{"mixmeta"}.}
  \item{\dots }{ further arguments passed to or from other methods.}
}

\details{
The model frame is produced by \code{\link{mixmeta}} when fitting the meta-analytical model, and stored in the \code{mixmeta} object if argument \code{model=TRUE}. Alternatively, the model frame is directly returned from a call to \code{\link{mixmeta}} with argument \code{method="model.frame"}. The method function \code{model.frame} simply extracts the saved model frame if available, or otherwise evaluates a call to \code{\link{mixmeta}} when \code{method="model.frame"}.

The method function \code{model.matrix} extracts the design matrix for the fixed-effects part of a fitted meta-analytical model. It first extract the model frame by calling \code{model.frame}, and then passes the call to the default method.

Note that the model frame of \code{mixmeta} models consist of terms for both the fixed and random-effects parts, the latter including also the grouping factors. This information can be used to reconstruct the proper model frame or matrix for each part.

These methods functions are similar to those provided for regression objects \code{\link{lm}} and \code{\link{lm}}.
}

\value{
For \code{model.frame}, a data.frame with special attributes (see the default method \code{\link{model.frame}}) and the additional class \code{"data.frame.mixmeta"}.

For \code{model.matrix}, the design matrix used to fit the model.
}

\author{Antonio Gasparrini <\email{antonio.gasparrini@lshtm.ac.uk}>}

\note{
The reason why these specific method functions are made available for class \code{mixmeta}, and in particular why a new class \code{"data.frame.mixmeta"} has been defined for model frames, lies in the special handling of missing values in multivariate meta-analysis models fitted with \code{\link{mixmeta}}. Methods \code{\link[=na.omit.data.frame.mixmeta]{na.omit}} and \code{\link[=na.exclude.data.frame.mixmeta]{na.exclude}} for class \code{"data.frame.mixmeta"} are useful for properly accounting for missing values when fitting these models.
}

\seealso{
See the default methods \code{\link{model.frame}} and \code{\link{model.matrix}}. See \code{\link[=na.omit.data.frame.mixmeta]{na.omit}} and \code{\link[=na.exclude.data.frame.mixmeta]{na.exclude}} on the handling of missing values. See \code{\link{mixmeta-package}} for an overview of the package and modelling framework.
}

\examples{
# RUN THE MODEL AND SUMMARIZE THE RESULTS
model <- mixmeta(cbind(PD,AL) ~ pubyear, S=berkey98[5:7], data=berkey98,
  method="ml")

# MODEL FRAME
model$model
model.frame(model)
update(model, method="model.frame")
class(model.frame(model))

# MODEL MATRIX
model.matrix(model)
}

\keyword{models}
\keyword{regression}
\keyword{multivariate}
\keyword{methods}
