\name{blup}
\alias{blup}

\title{ Best Linear Unbiased Predictions }

\description{
This is a generic function for generating best linear unbiased predictions (BLUPs) from the results of various fitting functions for meta-analytical models. The function invokes particular methods which depend on the \code{\link{class}} of the first argument. Currently, specific methods exist for several meta-analytical models in various packages: \code{\link[mixmeta]{blup.mixmeta}}, \code{\link[metafor]{blup.rma.uni}}, \code{\link[mvmeta]{blup.mvmeta}}, and \code{\link[dosresmeta]{blup.dosresmeta}}.
}

\usage{
blup(object, \dots)
}

\arguments{
  \item{object }{ a model object for which BLUPs are desired.}
  \item{\dots }{ further arguments passed to or from other methods.}
}

\details{
The generic method function \code{blup} calls specific method functions to produces (empirical) best linear unbiased predictions (BLUPs) from model objects.

These predictions are a shrunk version of unit-specific realizations, where unit-specific estimates borrow strength from the assumption of an underlying (potentially multivariate) distribution in a (usually hypothetical) population. The amount of shrinkage depends from the relative size of the within and between-unit covariance matrices.
}

\value{
The form of the value returned by \code{blup} depends on the class of its argument. See the documentation of the particular methods for details of what is produced by that method. Usually, the results consist of point estimates of BLUPs and optionally some measure of their uncertainty.
}

\references{
Verbeke G, Molenberghs G. \emph{Linear Mixed Models for Longitudinal Data}. Springer; 1997.

Sera F, Armstrong B, Blangiardo M, Gasparrini A (2019). An extended mixed-effects framework for meta-analysis.\emph{Statistics in Medicine}. 2019;38(29):5429-5444. [Freely available \href{http://www.ag-myresearch.com/2019_sera_statmed.html}{\bold{here}}].
}

\author{Antonio Gasparrini <\email{antonio.gasparrini@lshtm.ac.uk}> and Francesco Sera <\email{francesco.sera@lshtm.ac.uk}>}

\seealso{
Specific methods for various classes: \code{\link[mixmeta]{blup.mixmeta}}, \code{\link[metafor]{blup.rma.uni}}, \code{\link[mvmeta]{blup.mvmeta}}, and \code{\link[dosresmeta]{blup.dosresmeta}}.
}

\keyword{models}
\keyword{regression}
\keyword{multivariate}
\keyword{methods}
