\name{terms.mixmeta}
\alias{terms.mixmeta}


\title{ Extract Model Terms from mixmeta Objects }

\description{
These method function returns the terms object that defines meta-analytical models represented in objects of class \code{"mixmeta"}.
}

\usage{
\method{terms}{mixmeta}(x, type="fixed", \dots)
}

\arguments{
  \item{x }{ an object of class \code{"mixmeta"}.}
  \item{type }{ the type of terms. Either \code{"fixed"} or \code{"full"}. See Details.}
  \item{\dots }{ further arguments passed to or from other methods.}
}

\details{
The \code{terms} object is produced by \code{\link{mixmeta}} when fitting the meta-analytical model, and stored as an attribute of the \code{\link[=model.frame.mixmeta]{model.frame}}. Note that this object consists of terms for both the fixed and random-effects parts, the latter including also the grouping factors.

By using the default \code{type="fixed"}, this method function removes the random-effects terms. This can then be used, for instance, for creating the \code{\link[=model.matrix.mixmeta]{model.matrix}} for the fixed effects. Otherwise with \code{type="full"}, the full set of terms is returned.
}

\value{
An object of class \code{c("terms","formula")} which contains the terms representation of a symbolic meta-analytical model. See \code{\link{terms.object}} for its structure.
}

\author{Antonio Gasparrini <\email{antonio.gasparrini@lshtm.ac.uk}>}

\seealso{
See the  methods \code{\link[=model.frame.mixmeta]{model.frame}} and \code{\link[=model.matrix.mixmeta]{model.matrix}}.
}

\examples{
# RUN A MODEL
model <- mixmeta(effect, var, random= ~ 1|district/study, data=school)

# TERMS (FIXED AND FULL)
terms(model)
terms(model, "full")
attr(model.frame(model), "terms")
}

\keyword{models}
\keyword{regression}
\keyword{multivariate}
\keyword{methods}
