\name{predict.mixmeta}
\alias{predict.mixmeta}


\title{ Predicted Values from mixmeta Models }

\description{
This method function computes predictions from fitted univariate or multivariate meta-analytical models represented in objects of class \code{"mixmeta"}, optionally for a new set of predictor values in meta-regression models. Predictions are optionally accompanied by standard errors, confidence intervals or the entire (co)variance matrix of the predicted outcomes.
}

\usage{
\method{predict}{mixmeta}(object, newdata, se=FALSE, ci=FALSE, vcov=FALSE, ci.level=0.95,
  format, aggregate="stat", na.action=na.pass, \dots)
}

\arguments{
  \item{object }{ an object of class \code{"mixmeta"}.}
  \item{newdata }{ An optional data frame in which to look for variables values with which to predict from meta-regression models.}
  \item{se }{ logical switch indicating if standard errors must be included.}
  \item{ci }{ logical switch indicating if confidence intervals must be included.}
  \item{vcov }{ logical switch indicating if the (co)variance matrix must be included.}
  \item{ci.level }{ a numerical value between 0 and 1, specifying the confidence level for the computation of confidence intervals.}
  \item{format }{ the format for the returned results. See Value.}
  \item{aggregate }{ when \code{format="matrix"} and \code{se} or \code{ci} are required, the results may be aggregated by statistic or by outcome. See Value}
  \item{na.action }{ a function which indicates what should happen when the data contain \code{NA}s. The default to the value saved in \code{object}. See Note.}
  \item{\dots }{ further arguments passed to or from other methods.}

}

\details{
The method function \code{predict} produces predicted values from \code{mixmeta} objects, obtained by evaluating the original call to \code{\link{mixmeta}} in the frame \code{newdata}. For both fixed and random-effects models, estimated predictions are only based on the fixed part of the model, ignoring study-specific deviations, differently from \code{\link[=blup.mixmeta]{blup}}.

If \code{newdata} is omitted, the predictions are based on the data used for the fit. In that case how to handle predictions for units removed from estimation due to invalid missing pattern is determined by the \code{na.action} argument used in \code{\link{mixmeta}} to produce \code{object}. If \code{na.action=na.omit}, units excluded from estimation will not appear, whereas if \code{na.action=na.exclude} they will appear, with values set to \code{NA} for all the outcomes. This step is performed by \code{\link{napredict}}. See Notes.
}

\value{
The results may be aggregated in matrices (the default), or returned as lists, depending on the argument \code{format}. For multivariate models, the aggregation is ruled by the argument \code{aggregate}, and the results may be grouped by statistic or by outcome. If \code{vcov=TRUE}, lists are always returned.
}

\references{
Sera F, Armstrong B, Blangiardo M, Gasparrini A (2019). An extended mixed-effects framework for meta-analysis.\emph{Statistics in Medicine}. DOI: 10.1002/sim.8362 [Freely available \href{http://www.ag-myresearch.com/2019_sera_statmed.html}{\bold{here}}].
}

\author{Antonio Gasparrini <\email{antonio.gasparrini@lshtm.ac.uk}> and Francesco Sera <\email{francesco.sera@lshtm.ac.uk}>}

\note{
The definition of missing in model frames used for estimation in \code{\link{mixmeta}} is different than that commonly adopted in other regression models such as \code{\link{lm}} or \code{\link{glm}}. See info on \code{\link[=na.omit.data.frame.mixmeta]{missing values}} in \code{\link{mixmeta}}.
}

\seealso{
See \code{\link[=blup.mixmeta]{blup}} for best linear unbiased predictions. See the default method \code{\link{predict}}. See \code{\link{mixmeta-package}} for an overview of the package and modelling framework.
}

\examples{
# RUN THE MODEL
model <- mixmeta(cbind(PD,AL) ~ pubyear, S=berkey98[5:7], data=berkey98)

# PREDICTED FROM YEAR 1985 TO 1987, WITH LABELS
newdata <- data.frame(pubyear=1985:1987, row.names=1985:1987)

# AVERAGED OUTOCOMES AND SE
predict(model, newdata, se=TRUE)

# SAME AS ABOVE, AGGREGATED BY OUTCOME
predict(model, newdata, se=TRUE, aggregate="outcome")

# WITH VCOV, FORCED TO A LIST
predict(model, newdata, se=TRUE, vcov=TRUE, aggregate="outcome")
}

\keyword{models}
\keyword{regression}
\keyword{multivariate}
\keyword{methods}
