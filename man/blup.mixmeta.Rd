\name{blup.mixmeta}
\alias{blup.mixmeta}

\title{ Best Linear Unbiased Predictions from mixmeta Models }

\description{
This method function computes (empirical) best linear unbiased predictions from fitted random-effects meta-analytical models represented in objects of class \code{"mixmeta"}. Quantities can represent prediction of outcomes given both fixed and random effects, or just random-effects residuals from the fixed-effects estimates. Predictions are optionally accompanied by standard errors, prediction intervals or the entire (co)variance matrix of the predicted outcomes.
}

\usage{
\method{blup}{mixmeta}(object, se=FALSE, pi=FALSE, vcov=FALSE, pi.level=0.95, type="outcome",
  level, format, aggregate="stat", \dots)
}

\arguments{
  \item{object }{ an object of class \code{"mixmeta"}.}
  \item{se }{ logical switch indicating if standard errors must be included.}
  \item{pi }{ logical switch indicating if prediction intervals must be included.}
  \item{vcov }{ logical switch indicating if the (co)variance matrix must be included.}
  \item{pi.level }{ a numerical value between 0 and 1, specifying the confidence level for the computation of prediction intervals.}
  \item{type }{ the type of prediction. This can be either \code{outcome} (default) or \code{residual}. See Details.}
  \item{level }{ level of random-effects grouping for which predictions are to be computed. Default to the highest (inner) level, with 0 corresponding to fixed-effects predictions obtained through \code{\link[=predict.mixmeta]{predict}}.}
  \item{format }{ the format for the returned results. See Value.}
  \item{aggregate }{ when \code{format="matrix"} and \code{se} or \code{ci} are required, the results may be aggregated by statistic or by outcome. See Value.}
  \item{\dots }{ further arguments passed to or from other methods.}
}

\details{
The method function \code{blup} produces (empirical) best linear unbiased predictions from \code{mixmeta} objects. These can represent outcomes, given by the sum of fixed and random parts, or just random-effects residuals representing deviations from the fixed-effects estimated outcomes. In non-standard models with multiple hierarchies of random effects, the argument \code{level} can be used to determine the level of grouping for which predictions are to be computed.

These predictions are a shrunk version of unit-specific realizations, where unit-specific estimates borrow strength from the assumption of an underlying (potentially multivariate) distribution of outcomes or residuals in a (usually hypothetical) population. The amount of shrinkage depends from the relative size of the within and between-unit covariance matrices reported as components \code{S} and \code{Psi} in \code{mixmeta} objects (see \code{\link{mixmetaObject}}).

Fixed-effects models do not assume random effects, and the results of \code{blup} for these models are identical to \code{\link[=predict.mixmeta]{predict}} (for \code{type="oucome"}) or just 0's (for \code{type="residuals"}).

How to handle predictions for units removed from estimation due to invalid missing pattern is determined by the \code{na.action} argument used in \code{\link{mixmeta}} to produce \code{object}. If \code{na.action=na.omit}, units excluded from estimation will not appear, whereas if \code{na.action=na.exclude} they will appear, with values set to \code{NA} for all the outcomes. This step is performed by \code{\link{napredict}}. See Note below.

In the presence of missing values in the outcomes \code{y} of the fitted model, correspondent values of point estimates and covariance terms are set to 0, while the variance terms are set to \code{1e+10}. In this case, in practice, the unit-specific estimates do not provide any information (their weight is virtually 0), and the prediction tends to the value returned by \code{\link[=predict.mixmeta]{predict}} with \code{interval="prediction"}, when applied to a new but identical set of predictors. See also Note below.
}

\value{
(Empirical) best linear unbiased predictions of outcomes or random-effects residuals. The results may be aggregated in matrices (the default), or returned as lists, depending on the argument \code{format}. For multivariate models, the aggregation is ruled by the argument \code{aggregate}, and the results may be grouped by statistic or by outcome. If \code{vcov=TRUE}, lists are always returned.
}

\references{
Sera F, Armstrong B, Blangiardo M, Gasparrini A (2019). An extended mixed-effects framework for meta-analysis.\emph{Statistics in Medicine}. 2019;38(29):5429-5444. [Freely available \href{http://www.ag-myresearch.com/2019_sera_statmed.html}{\bold{here}}].

Verbeke G, Molenberghs G. \emph{Linear Mixed Models for Longitudinal Data}. Springer; 1997.
}

\author{Antonio Gasparrini <\email{antonio.gasparrini@lshtm.ac.uk}> and Francesco Sera <\email{francesco.sera@lshtm.ac.uk}>}

\note{
The definition of missing in model frames used for estimation in \code{\link{mixmeta}} is different than that commonly adopted in other regression models such as \code{\link{lm}} or \code{\link{glm}}. See info on \code{\link[=na.omit.data.frame.mixmeta]{missing values}} in \code{\link{mixmeta}}.

Differently from \code{\link[=predict.mixmeta]{predict}}, this method function computes the predicted values in the presence of partially missing outcomes. Interestingly, BLUPs for missing outcomes may be slightly different than predictions returned by \code{\link[=predict.mixmeta]{predict}} on a new but identical set of predictors, as the BLUP also depends on the random part of the model. Specifically, the function uses information from the random-effects (co)variance to predict missing outcomes given the observed ones.
}

\seealso{
See \code{\link[=predict.mixmeta]{predict}} for standard predictions. See \code{\link{mixmeta-package}} for an overview of the package and modelling framework.
}

\examples{
# RUN THE MODEL
model <- mixmeta(cbind(PD,AL) ~ 1, S=berkey98[5:7], data=berkey98)

# ONLY BLUP
blup(model)

# BLUP AND SE
blup(model, se=TRUE)

# SAME AS ABOVE, AGGREGATED BY OUTCOME, WITH PREDICTION INTERVALS
blup(model, se=TRUE, pi=TRUE, aggregate="outcome")

# WITH VCOV, FORCED TO A LIST
blup(model, se=TRUE, pi=TRUE, vcov=TRUE, aggregate="outcome")

# PREDICTING ONLY THE RANDOM-EFFECT RESIDUALS
blup(model, type="residual")
}

\keyword{models}
\keyword{regression}
\keyword{multivariate}
\keyword{methods}
