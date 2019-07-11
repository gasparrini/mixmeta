\name{qtest}
\alias{qtest}


\title{ Cochran Q Test of Heterogeneity }

\description{
This is a generic function to perform a Cochran Q test of (residual) heterogeneity. The function invokes particular \code{\link{methods}} which depend on the \code{\link{class}} of the first argument. Currently, specific methods exist for several meta-analytical models in various packages: \code{\link[mixmeta]{qtest.mixmeta}}, \code{\link[mvmeta]{qtest.mvmeta}}, and \code{\link[dosresmeta]{qtest.dosresmeta}}.
}

\usage{
qtest(object, \dots)
}

\arguments{

  \item{object }{ an object for which the test is desired}
  \item{\dots}{ further arguments passed to specific methods.}
}

\details{
The test assesses the null hypothesis that the variability in the distribution of the outcomes is explained only in terms of within-unit estimation errors. This corresponds to a test on the hypotesis that there is no variation attributable to random-effects terms.
}

\value{
Returned values depend on the specific class. Usually, the results of the test.
}

\references{
Cochran WG (1950). The comparison of percentages in matched samples". \emph{Biometrika}. \bold{37}(3/4):256--266.

Sera F, Gasparrini A (2019). An extended mixed-effects framework for meta-analysis.\emph{Statistics in Medicine}. In submission.
}

\author{Antonio Gasparrini <\email{antonio.gasparrini@lshtm.ac.uk}>}

\seealso{
Specific methods for various classes: \code{\link[mixmeta]{qtest.mixmeta}}, \code{\link[mvmeta]{qtest.mvmeta}}, and \code{\link[dosresmeta]{qtest.dosresmeta}}.}

\keyword{methods}
\keyword{htest}
