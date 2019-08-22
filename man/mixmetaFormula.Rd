\name{mixmetaFormula}
\alias{mixmetaFormula}


\title{ Formulae in mixmeta Models }

\description{
Meta-analytical models fitted with \code{\link{mixmeta}} are defined by specific formulae in its arguments \code{formula} and \code{random}. The formulae offer compact symbolic expressions with form \code{y ~ x + z}, where the response \code{y} in the left-hand side of the operator \code{~} is modelled in terms of meta-predictors \code{x} and \code{z} in the right-hand side. Terms are separated by \code{+}, and additional syntactic operators and existing functions can be used within a formula to specify transformations such as categorization and interactions, among others, as in standard formula expressions (see \code{\link[stats]{formula}} for details). The usage of formulae in \pkg{mixmeta} for the random-effects part follows closely the definition in the the \pkg{nlme} package.
}

\section{Fixed-effects formula}{
The argument \code{formula} of \code{\link{mixmeta}} defines the fixed-effects part. Models for meta-analysis with no meta-predictors can be specified using the form \code{y ~ 1}, or alternatively including only the term \code{y} (in this case, the formula is reconstructed internally). Multivariate models can be defined by using a matrix-type \code{y}, with columns as multiple outcomes, or directly in the formula with form \code{cbind(y1 + y2) ~ 1}. In  meta-regression models, other terms are added in the right-hand side of the formula as \code{y ~ x1 + ... + xp}, defining the linear meta-predictor. In multivariate meta-regression, the same linear predictor is specified for each outcome.
}

\section{Random-effects formula or formulae}{
The argument \code{random} of \code{\link{mixmeta}} defines the random-effects part. When this is not specified, it is assumed that each row of data is from an independent study and assigned to a different group, as in standard meta-analytical models. If provided, this is usually represented by a one-sided formula whose basic random-intercept form is \code{~ 1 | g}. The term \code{g} at the right-hand side of the special operator \code{|} is the grouping factor, always required in a single random-effects formula. A more complex random-effects part can be also specified by \code{~ z1 + ... + zq | g}, where the terms in the left-hand side defines a linear meta-predictor, with syntax identical to the usual formulae.

The argument \code{random} also accepts a list of one-sided formulae referring to multiple random-effects levels in multilevel meta-analytical models. In this case, levels are assumed to be nested in the order of the list, from the lowest (outer) to the highest (inner), consistently with the grouping factors. These are usually defined by different terms in the right-hand side of the \code{|} operator, although in the list form can also be provided as names of the list components. This information is used internally to reconstruct the grouping structure and the random-effects design matrices. Each level can have different linear predictors, but if these are identical across levels the random-effects part can be defined by a single equation \code{~ z | g1 / g2}, where the special operator \code{/} separates the grouping factors \code{g2} nested in \code{g1}.
}

\references{
Sera F, Gasparrini A (2019). An extended mixed-effects framework for meta-analysis.\emph{Statistics in Medicine}. DOI: 10.1002/sim.8362 [Freely available \href{http://www.ag-myresearch.com/2019_sera_statmed.html}{\bold{here}}].

}

\author{Antonio Gasparrini <\email{antonio.gasparrini@lshtm.ac.uk}> and Francesco Sera <\email{francesco.sera@lshtm.ac.uk}>}

\seealso{
See \code{\link{mixmeta}}. See \code{\link{formula}} for standard regression formulae. See \code{\link{mixmeta-package}} for an overview of this modelling framework.
}

\examples{
# STANDARD RANDOM-EFFECTS META-ANALYSIS (WITH DIFFERENT SYNTAXES)
mixmeta(logor, logorvar, data=bcg)
mixmeta(logor ~ 1, logorvar, data=bcg)

# META-REGRESSION
mixmeta(logor ~ ablat, logorvar, data=bcg)

# MULTIVARIATE MODEL
model <- mixmeta(cbind(PD,AL) ~ pubyear, S=berkey98[5:7], data=berkey98)

# NON-STANDARD MODEL: REPEATED MEASURED WITHING THE SAME GROUPS
mixmeta(effect, var, random= ~ 1|district, data=school)
mixmeta(absrisk, var, random= ~ 1|trial, data=thrombolytic)

# NON-STANDARD MODEL: MORE COMPLEX RANDOM-EFFECTS PREDICTOR
mixmeta(logOR~time, logORvar, random= ~ I(time-15)|study, data=gliomas)

# MULTILEVEL MODEL (WITH DIFFERENT SYNTAXES)
mixmeta(effect, var, random= ~ 1|district/study, data=school)
mixmeta(effect, var, random=list(~ 1|district, ~ 1|study), data=school)
mixmeta(effect, var, random=list(district = ~ 1, study = ~ 1), data=school)
}

\keyword{models}
\keyword{regression}
\keyword{multivariate}

