\name{na.omit.data.frame.mixmeta}
\alias{na.omit.data.frame.mixmeta}
\alias{na.exclude.data.frame.mixmeta}


\title{ Handling Missing Values in mixmeta Models }

\description{
These method functions exclude rows corresponding to units with invalid missing pattern from model frames of class \code{"data.frame.mixmeta"}. This guarantees the correct handling of missing values while fitting meta-analytical models.
}

\usage{
\method{na.omit}{data.frame.mixmeta}(object, \dots)

\method{na.exclude}{data.frame.mixmeta}(object, \dots)
}

\arguments{
  \item{object }{ an object of class \code{"data.frame.mixmeta"}.}
  \item{\dots }{ further arguments passed to or from other methods.}
}

\details{
A model frame of class \code{"data.frame.mixmeta"} is produced by \code{\link{mixmeta}}. A call to \code{na.omit} or \code{na.exclude} removes from the model frame the rows corresponding to studies with invalid missing pattern. In addition, a \code{na.action} attribute is added to the model frame, namely a numeric vector corresponding to the removed rows and class \code{"omit"} or \code{"exclude"}, respectively. This information is used by \code{\link{naresid}} and \code{\link{napredict}} to deal with missing values in functions such as \code{\link{fitted}}, \code{\link{residuals}}, \code{\link[=predict.mixmeta]{predict}} and \code{\link[=blup.mixmeta]{blup}}, among others.

The definition of missing, identifying an invalid missing pattern, is different in meta-analytical models performed through \code{\link{mixmeta}} if compared to other regression functions such as \code{\link{lm}} or \code{\link{glm}}, in particular for the multivariate case. Specifically, while a unit is removed if at least an observation for one predictor is missing, partially missing outcomes do not prevent the unit to contribute to estimation (see \code{\link{mixmeta}}). Specific methods \code{na.omit} and \code{na.exclude} for class \code{"data.frame.mixmeta"} allow  this different definition.
}

\value{
These functions returns the model frame \code{object} with rows corresponding to units with invalid missing pattern being removed. They also add the related \code{na.action} attribute as explained above.
}

\author{Antonio Gasparrini <\email{antonio.gasparrini@lshtm.ac.uk}>}

\seealso{
See \code{\link{na.action}}, \code{\link{naresid}} and \code{\link{napredict}}. See \code{\link[=model.frame.mixmeta]{model.frame}}.
See \code{\link{mixmeta-package}} for an overview of the package and modelling framework.
}

\examples{
# INPUT MISSING VALUES IN PREDICTOR AND ONE RESPONSE
data <- berkey98
data[2,1] <- data[4,3] <- NA
data

# RUN THE MODEL
model <- mixmeta(cbind(PD,AL) ~ pubyear, S=data[5:7], data=data, method="ml")

# SUMMARIZE: NOTE THE NUMBER OF STUDIES AND OBSERVATIONS
summary(model)
df.residual(model)

# EXTRACT THE MODEL FRAME WITH na.pass
model.frame(model, na.action="na.pass")
# EXTRACT THE MODEL FRAME WITH na.omit (DEFAULT)
model.frame(model, na.action="na.omit")

# COMPARE WITH DEFAULT METHOD FOR na.omit
frame <- model.frame(model, na.action="na.pass")
na.omit(frame)
class(frame)
class(frame) <- "data.frame"
na.omit(frame)

# WITH na.exclude
residuals(model)
residuals(update(model, na.action="na.exclude"))
}

\keyword{models}
\keyword{regression}
\keyword{multivariate}
\keyword{methods}
\keyword{manip}
