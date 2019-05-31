\name{berkey98}
\alias{berkey98}
\docType{data}

\title{ Five Published Trials on Periodontal Disease }

\description{
The dataset contains the results of 5 published trials comparing surgical and non-surgical treatments for medium-severity periodontal disease, one year after treatment. The 2 estimated outcomes are average improvements (surgical minus non-surgical, in mm) in probing depth (PD) and attachment level (AL).
}

\usage{berkey98}

\format{
  A data frame with 5 observations on the following 7 variables:
  \describe{
    \item{\code{pubyear}}{ publication year of the trial.}
    \item{\code{npat}}{ number of patients included in the trial.}
    \item{\code{PD}}{ estimated improvement of surgical versus non-surgical treatments in probing depth (mm).}
    \item{\code{AL}}{ estimated improvement of surgical versus non-surgical treatments in attachment level (mm).}
    \item{\code{var_PD}}{ variance of the estimated outcome for \code{PD}.}
    \item{\code{cov_PD_AL}}{ covariance of the estimated outcomes for \code{PD} and \code{AL}.}
    \item{\code{var_AL}}{ variance of the estimated outcome for \code{AL}.}
  }
  Row names specify the author of the paper reporting the results of each trial.
}

\source{
Berkey CS, et al. (1998). Meta-analysis of multiple outcomes by regression with random effects. \emph{Statistics in Medicine}. \bold{17}:2537--2550.

Berkey CS., et al. (1995). Multiple-outcomes meta-analysis of treatments for periodontal disease. \emph{Journal of Dental Research}. \bold{74}(4):1030--1039.
}

\examples{
### REPRODUCE THE RESULTS IN BERKEY ET AL. (1998)

# INSPECT THE DATA
berkey98

# FIXED-EFFECTS
year <- berkey98$pubyear - 1983
mod1 <- mixmeta(cbind(PD,AL) ~ year, S=berkey98[5:7], data=berkey98,
  method="fixed")
print(summary(mod1), digits=3)

# GLS MODEL (VARIANCE COMPONENTS)
mod2 <- mixmeta(cbind(PD,AL) ~ year, S=berkey98[5:7], data=berkey98,
  method="vc",  control=list(vc.adj=FALSE))
print(summary(mod2), digits=3)
round(mod2$Psi, 3)

# ML MODEL
mod3 <- mixmeta(cbind(PD,AL) ~ year, S=berkey98[5:7], data=berkey98, method="ml")
print(summary(mod3), digits=3)
round(mod3$Psi, 3)
}

\keyword{datasets}
