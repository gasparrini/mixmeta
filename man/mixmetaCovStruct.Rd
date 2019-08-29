\name{mixmetaCovStruct}
\alias{mixmetaCovStruct}


\title{ (Co)variance Structures for mixmeta Models }

\description{
Alternative options for the (co)variance structure of the random effects random effects in meta-analytical models, usually defined through the argument \code{bscov} of the function \code{\link{mixmeta}}.
}

\section{Options}{
Assuming a meta-analysis or meta-regression based on \eqn{k} outcomes, for each grouping level with \eqn{q} random-effects predictors the matrix can be specified in various forms listed below. For multivariate models with multiple predictors, the order implies a sequence of \eqn{q} parameters for each \eqn{k} outcomes. These are the options:

  \itemize{
    \item \code{unstr}: an unstructured form for a general positive-definite matrix. The matrix is represented by \eqn{kq(kq+1)/2} unrestricted parameters defined as the upper triangular entries of its Cholesky decomposition.
    \item \code{diag}: a diagonal positive-definite matrix. The matrix is represented by \eqn{kq} unrestricted parameters defined as the logarithm of the diagonal values.
    \item \code{id}: a multiple of the identity positive-definite matrix. The matrix is represented by a single unrestricted parameter defined as the logarithm of the diagonal value.
    \item \code{cs}: a positive-definite matrix with compound symmetry structure. The matrix is represented by 2 unrestricted parameters defined as the logarithm of the identical diagonal value and the transformed correlation. The latter is parameterized so to obtain a correlation value between \eqn{-1/(kq-1)} and 1, in order to ensure positive-definiteness.
    \item \code{hcs}: a positive-definite matrix with heterogeneous compound symmetry structure. The matrix is represented by \eqn{kq+1} unrestricted parameters defined as the logarithm of the diagonal values and the transformed correlation. The latter is parameterized so to obtain a correlation value between \eqn{-1/(kq-1)} and 1, in order to ensure positive-definiteness.
    \item \code{ar1}: a positive-definite matrix with autoregressive structure of first order. The matrix is represented by 2 unrestricted parameters defined as the logarithm of the identical diagonal value and the logistic transformed correlation. The latter is parameterized so to obtain a correlation value between -1 and 1.
   \item \code{har1}: a positive-definite matrix with heterogeneous autoregressive structure of first order. The matrix is represented by \eqn{kq+1} unrestricted parameters defined as the logarithm of the diagonal value and the logistic transformed correlation. The latter is parameterized so to obtain a correlation value between -1 and 1.
    \item \code{prop}: a positive-definite matrix proportional to that provided by the user through the argument \code{Psifix} in the control list (see \code{\link{mixmeta.control}}). The matrix is represented by 1 unrestricted parameter defined as the logarithm of the multiplier.
    \item \code{cor}: a positive-definite matrix with correlation structure provided by the user through the argument \code{Psifix} (with \code{\link{cov2cor}}) in the control list (see \code{\link{mixmeta.control}}). The matrix is represented by \eqn{k} unrestricted parameters defined as the logarithm of the diagonal values.
    \item \code{fixed}: a known matrix provided by the user through the argument \code{Psifix} in the control list (see \code{\link{mixmeta.control}}). The matrix is known and no parameters are needed to represent it.
  }
}

\section{Details}{
Structures other than \code{unstr} are only available for models estimated through (restricted) maximum likelihood.

The unrestricted parameters defining the random-effects (co)variance matrix (or matrices for multilvel models) are estimated in the iterative optimization algorithm (see \code{\link{mixmeta.ml}}). Although rarely needed and not recommeded, the user can provided a starting value of the (co)variance matrix, from which the parameters are derived (see \code{\link{mixmeta.control}}).
}

\references{
Sera F, Armstrong B, Blangiardo M, Gasparrini A (2019). An extended mixed-effects framework for meta-analysis.\emph{Statistics in Medicine}. DOI: 10.1002/sim.8362 [Freely available \href{http://www.ag-myresearch.com/2019_sera_statmed.html}{\bold{here}}].

Pinheiro JC and Bates DM (2000). \emph{Mixed-Effects Models in S and S-PLUS}. New York, Springer Verlag.
}

\author{Antonio Gasparrini <\email{antonio.gasparrini@lshtm.ac.uk}> and Francesco Sera <\email{francesco.sera@lshtm.ac.uk}>}

\section{Note}{
The choice of structures can affect the performance of the optimization procedure, determining forms of likelihood surfaces which induce convergence to local maxima. In particular, structures such as multiple of identity or proportional to a fixed matrix are based on strong assumptions and should be used with caution.
}

\seealso{
See \code{\link{mixmeta}}. See \code{\link{lm}} or \code{\link{glm}} for standard regression functions. See \code{\link{mixmeta-package}} for an overview of this modelling framework.
}

\examples{
# UNSTRUCTURED AND STRUCTURED BETWEEN-STUDY (CO)VARIANCE
y <- as.matrix(smoking[11:13])
S <- as.matrix(smoking[14:19])
mod1 <- mixmeta(y, S)
summary(mod1)
mod1$Psi

# DIAGONAL
mod2 <- mixmeta(y, S, bscov="diag")
summary(mod2)
mod2$Psi

# HETEROGENEOUS COMPOUND SYMMETRY
mod3 <- mixmeta(y, S, bscov="hcs")
summary(mod3)
mod3$Psi

# PROPORTIONAL
mod4 <- mixmeta(y, S, bscov="prop", control=list(Psifix=diag(3)+1))
summary(mod4)
mod4$Psi

# CORRELATION
Psicor <- matrix(0.2, 3, 3) ; diag(Psicor) <- 1
mod5 <- mixmeta(y, S, bscov="cor", control=list(Psifix=Psicor))
summary(mod5)
mod5$Psi
}

\keyword{models}
\keyword{regression}
\keyword{multivariate}

