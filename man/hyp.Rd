\name{hyp}
\alias{hyp}
\docType{data}

\title{ Ten Studies Assessing an Hypertension Treatment }

\description{
The dataset contains the results of ten studies that assess the effectiveness of hypertension treatment for lowering blood pressure. Each study provides complete data on two treatment effects, the difference in systolic blood pressure (SBP) and diastolic blood pressure (DBP) between the treatment and the control groups, where these differences are adjusted for the participants' baseline blood pressures. The within-study correlations of the two outcomes are known. Some trials are conducted on patients with isolated systolic hypertension (ISH).
}

\usage{hyp}

\format{
  A data frame with 10 observations on the following 7 variables:
  \itemize{
    \item{\code{study}: study ID.}
    \item{\code{sbp, sbp_se}: estimated difference and its standard error in systolic blood pressure.}
    \item{\code{dbp, dbp_se}: estimated difference and its standard error in diastolic blood pressure.}
    \item{\code{rho}: within-study correlation between the estimated differences in systolic and diastolic blood pressure.}
    \item{\code{ish}: indicator for studies on patients with isolated systolic hypertension.}
  }
}

\source{
Jackson D, White IR, Riley RD (2013). A matrix based method of moments for fitting the multivariate random effects model for meta-analysis and meta-regression. \emph{Biometrical Journal}. \bold{55}(2):231--45.

Sera F, Gasparrini A (2019). An extended mixed-effects framework for meta-analysis.\emph{Statistics in Medicine}. In submission.
}

\note{
The standard errors for the two outcomes are wrongly reported as variances in the original article by Jackson and colleagues (2013).
}

\examples{
### REPRODUCE THE RESULTS IN JACKSON ET AL. (2013)

# INSPECT THE DATA
hyp

# INPUT THE CORRELATION (CAN ALSO BE INPUTTED DIRECTLY, SEE BELOW)
(S <- inputcov(hyp[c("sbp_se", "dbp_se")], cor=hyp$rho))
# CHECK WITH THE FIRST STUDY
cov2cor(xpndMat(S[1,]))

# META-ANALYSIS, REML MODEL
mod1 <- mixmeta(cbind(sbp,dbp), S=S, data=hyp)
print(summary(mod1), digits=2)
round(mod1$Psi,2)

# META-ANALYSIS, REML MODEL (INPUTTING THE CORRELATION DIRECTLY)
mod2 <- mixmeta(cbind(sbp,dbp), S=cbind(sbp_se,dbp_se)^2, data=hyp,
  control=list(Scor=hyp$rho))
print(summary(mod2), digits=2)

# META-ANALYSIS, MM MODEL
mod3 <- mixmeta(cbind(sbp,dbp), S=S, data=hyp, method="mm")
print(summary(mod3), digits=2)
round(mod3$Psi,2)

# META-REGRESSION, REML MODEL
mod4 <- mixmeta(cbind(sbp,dbp) ~ ish, S=S, data=hyp)
print(summary(mod4), digits=2)

# META-REGRESSION, MM MODEL
mod5 <- mixmeta(cbind(sbp,dbp) ~ ish, S=S, data=hyp, method="mm")
print(summary(mod5), digits=2)
}

\keyword{datasets}
