\name{mixmeta.ml}
\alias{mixmeta.ml}
\alias{mixmeta.reml}


\title{ ML and REML Estimators for mixmeta Models }

\description{
These functions implement maximum likelihood (ML) and restricted maximum likelihood (REML) estimators for random-effects meta-analysis and meta-regression, including standard univariate models, and non-standard extensions such as multivariate, multilevel, longitudinal, and dose-response models. These functions are meant to be used internally and not directly run by the users.
}

\usage{
mixmeta.ml(Xlist, Zlist, ylist, Slist, nalist, rep, k, q, nall, bscov, control, ...)

mixmeta.reml(Xlist, Zlist, ylist, Slist, nalist, rep, k, q, nall, bscov, control, ...)
}

\arguments{
Assuming a meta-analysis or meta-regression based on \eqn{n} units aggregated within \eqn{m} (outer-level) groups, \eqn{k} outcomes, \eqn{p} fixed-effects predictors, and \eqn{q} random-effects predictors:

  \item{Xlist }{ a \eqn{m}-dimensional list of group-specific design matrices for the fixed-effects part of the model. Rows corresponding to missing outcomes have been excluded.}
  \item{Zlist }{ a \eqn{m}-dimensional list of group-specific design matrices for the random-effects part of the model. Each element of this list represents a list of matrices correponding to the (optionally multiple) grouping levels of random effects. In each matrix, rows corresponding to missing outcomes have been excluded.}
  \item{ylist }{ a \eqn{m}-dimensional list of group-specific vectors of estimated outcomes. Entries corresponding to missing outcomes have been excluded.}
  \item{Slist }{ a \eqn{m}-dimensional list of within-group (co)variance matrices of estimated outcomes. Rows and columns corresponding to missing outcomes have been excluded.}
  \item{nalist }{ a \eqn{m}-dimensional list of group-specific logical vectors, identifying missing outcomes.}
  \item{rep }{ matrix with \eqn{m} rows where each column identifies the number of repetitions (number of groups) for each grouping level. The first column (outer level) is by definition a vector of 1's.}
  \item{k, q, nall }{ number of outcomes, number of random-effects predictors (including the intercept), total number of observations (excluding missing), respectively. While usually all are scalars, in the case of multilevel models \code{q} can be a numeric vector representing the number of predictors for each level.}
  \item{bscov }{ a character vector defining the structure of the (co)variance matrix for each level or random effects. See \code{\link{mixmeta}}.}
  \item{control }{ list of parameters for controlling the fitting process, usually internally set to default values by \code{\link{mixmeta.control}}.}
  \item{\dots }{ further arguments passed to or from other methods. Currently not used.}
}

\details{
The estimation involves \eqn{kp} fixed-effects coefficients and random-effects parameters, whose number depends on the number of grouping levels and, for each of them, on the chosen \code{\link[=mixmetaCovStruct]{(co)variance structure}} for the between-study (co)variance matrices. A maximum of \eqn{kq(kq+1)/2} parameters are needed in the case of or single-level models with unstructured form for the random-effects (co)variance matrix.

(Restricted) maximum likelihood estimators implemented in \pkg{mixmeta} rely on two iterative algorithms: \code{\link[=ml.igls]{(R)IGLS}} and \code{\link[=ml.newton]{quasi-Newton}} iterative methods. The former implements a (restricted) iterative generalized least squares method, while the latter is based on a Newton-type maximization routine using specific \code{\link[=ml.loglik.fn]{likelihood functions}}. The default estimation method is based on a hybrid procedure, with few runs of of the (R)IGLS algorithm and then quasi-Newton iterations until convergence. This approach is optimal in exploiting the properties of both algorithms, with (R)IGLS being robust to the choice of initial values and quick in getting near the maximum, while the quasi-Newton is fast to converge from that point. Full (R)IGLS or quasi-Newton methods can be alternatively selected using the \code{control} argument of \code{\link{mixmeta}} (see \code{\link{mixmeta.control}}). Follow the links above for details on each iterative algorithm.

Both estimation algorithms adopt a profiled (or concentrated) approach, where the optimization is expressed only in terms of the random-effects parameters. Cholesky and and QR decompositions are used for computational stability and efficiency, and for assuring the positive-definiteness of the estimated between-study (co)variance matrix. The method is described in details in Gasparrini and collaborators (2012) (see references below).
}

\value{
These functions return an intermediate list object, with some components then processed and some others added later within \code{\link{mixmeta.fit}} and \code{\link{mixmeta}} to finalize an object of class \code{"mixmeta"}. See \code{\link{mixmetaObject}}.
}

\references{
Gasparrini A, Armstrong B, Kenward MG (2012). Multivariate meta-analysis for non-linear and other multi-parameter associations. \emph{Statistics in Medicine}. \bold{31}(29):3821--3839. [Freely available \href{http://www.ag-myresearch.com/2012_gasparrini_statmed.html}{\bold{here}}].

Pinheiro JC and Bates DM (2000). \emph{Mixed-Effects Models in S and S-PLUS}. New York, Springer Verlag.

Lindstrom MJ and Bates DM (1988). Newton-Raphson and EM algorithms for linear mixed-effects models for repeated-measures data. \emph{Journal of the American Statistical Association}. \bold{83}(404):1014--1022.

Goldstein H (1986). Multilevel mixed linear model analysis using iterative generalized least squares. \emph{Biometrika}. \bold{73}(1):43.

Goldstein H (1992). Efficient computational procedures for the estimation of parameters in multilevel models based on iterative generalized least squares. \emph{Computational Statistics \& Data Analysis}. \bold{13}(1):63--71.
}

\author{Antonio Gasparrini <\email{antonio.gasparrini@lshtm.ac.uk}> and Francesco Sera <\email{francesco.sera@lshtm.ac.uk}>}

\note{
As stated earlier, these functions are called internally by \code{\link{mixmeta.fit}}, and are not meant to be used directly. In particular, their code does not contain any check on the arguments provided, which are expected in specific formats. The functions are not exported in the namespace, and only documented for completeness.

The arguments above are prepared by \code{\link{mixmeta.fit}} from its arguments \code{X}, \code{Z}, \code{y}, \code{S}, \code{groups}, and \code{bscov}. The list structure, although requiring more elaborate coding, is computationally more efficient, as it avoids the specification of sparse block-diagonal matrices, especially for meta-analysis involving a large number of studies.

Some parameters of the fitting procedures are determined by the \code{control} argument, with default set by \code{\link{mixmeta.control}}. No missing values are accepted in the fitting functions. See details on \code{\link[=na.omit.data.frame.mixmeta]{missing values}}.
}

\seealso{
See \code{\link{mixmeta}} for the general usage of the functions. See \code{\link{mixmeta.control}} to determine specific parameters of the fitting procedures. Use the triple colon operator ('\code{:::}') to access the code of the internal functions, such as \code{glsfit}. See \code{\link{mixmeta-package}} for an overview of the package and modelling framework.
}

\examples{
# REML ESTIMATOR: UNIVARIATE MODEL
mod1 <- mixmeta(yC, S=SCC, data=smoking)
summary(mod1)

# ML ESTIMATOR: MULTIVARIATE MODEL
year <- berkey98$pubyear - 1983
mod2 <- mixmeta(cbind(PD,AL) ~ year, S=berkey98[5:7], data=berkey98,method="ml")
print(summary(mod2), digits=3)
round(mod2$Psi,3)

# STRUCTURED BETWEEN-STUDY (CO)VARIANCE
y <- as.matrix(fibrinogen[2:5])
S <- as.matrix(fibrinogen[6:15])
mod3 <- mixmeta(y, S, bscov="hcs")
summary(mod3)

# MULTILEVEL MODEL
mod4 <- mixmeta(effect, var, random= ~ 1|district/study, data=school)
summary(mod4)

# LONGITUDINAL MODEL
mod5 <- mixmeta(logOR~time, S=logORvar, random=~I(time-15)|study, bscov="diag",
  method="ml", data=gliomas)
summary(mod5)
}


\keyword{models}
\keyword{regression}
\keyword{multivariate}
