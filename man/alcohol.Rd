\name{alcohol}
\alias{alcohol}
\docType{data}

\title{ Alcohol Intake and Colorectal Cancer }

\description{
The dataset contains the data on 8 cohort studies participating in the Pooling Project of Prospective Studies of Diet and Cancer. A total of 3,646 cases and 2,511,424 person-years were included in the analysis. Each study estimated the incidence relative rate in different categories of alcohol intake while controlling for a set of potential confounders, using non-drinkers as the reference. The categories where then converted in a dose by assigning to each the median value of individual consumptions, with studies reporting estimates at different levels in a continuous scale.
}

\usage{alcohol}

\format{
  A data frame with 48 observations on the following 7 variables:
  \itemize{
    \item{\code{id}: label for each study, derived from the first author's name.}
    \item{\code{type}: code for study design (cohort estimating incidence rate).}
    \item{\code{dose}: assigned dose level (gr/day of alcohol intake).}
    \item{\code{cases}: number of cases for each dose category.}
    \item{\code{peryears}: amount of person-time for each dose category.}
    \item{\code{logrr}: estimated logarithm of the incidence relative rate.}
    \item{\code{se}: standard error of the estimates.}
  }
}

\details{
The data are stored in a \emph{long} format, with each record reporting the information for each dose categories and studies including multiple records. The reference category for each study included, although the log-RR is fixed to 0 with no standard error (comparing the category with itself). The information on these reference categories is needed to compute the approximate correlations between estimates in the same study.
}

\source{
Crippa A, et al (2018). One-stage dose-response meta-analysis for aggregated data. \emph{Statistical Methods in Medical Research}. DOI: 10.1177/0962280218773122.

Crippa A, Orsini N (2016). Multivariate dose-response meta-analysis: The dosresmeta R package. \emph{Journal of Statistical Software}. \bold{72}(1):1--15.

Orsini N, et al (2012). Meta-analysis for linear and nonlinear dose-response relations: examples, an evaluation of approximations, and software. \emph{American Journal of Epidemiology}. \bold{175}(1):66--73.
}

\note{
The data provide an example of application of dose-response meta-analysis, with repeated measurements of the effect size associated to different doses within each study. This requires a modelling structure that accounts for both within and between-study correlations of repeated measurements. The within-study correlations are usually reconstructed from published data using specific methods. Results can be compared with those reported by Crippa and Orsini (2016) and Orsini and colleagues (2012), although they are not identical: while the original analysis used a two-stage approach, the modelling framework applied here follows the more recent one-stage dose-response meta-analysis proposed by Crippa and collegues (2018).

The dataset is also available in the same format in the dataframe \code{alcohol_crc} of the package \pkg{dosresmeta}.
}

\examples{
### REPRODUCE THE RESULTS IN CRIPPA ET AL (2016) AND ORSINI ET AL (2012)

# LOAD THE PACKAGE dosresmeta AND splines
library(dosresmeta) ; library(splines)

# COMPUTE THE WITHIN-STUDY CORRELATIONS EXCLUDING THE REFERENCE
addS <- lapply(split(alcohol, alcohol$id), function(x)
  covar.logrr(y=logrr, v=se^2, cases=cases, n=peryears, type=type, data=x))
sub <- subset(alcohol, !is.na(se))

# NOT ACCOUNTING FOR WITHIN-STUDY CORRELATIONS
nocor <- mixmeta(logrr ~ 0 + dose, S=se^2, random= ~ 0 + dose|id, data=sub,
  method="ml")
summary(nocor)

# ACCOUNTING FOR WITHIN-STUDY CORRELATIONS
lin <- mixmeta(logrr ~ 0 + dose, random= ~ 0 + dose|id, data=sub, method="ml",
  control=list(addSlist=addS))
summary(lin)

# NATURAL SPLINES
knots <- quantile(alcohol$dose, c(5,35,65,95)/100)
spl <- ns(sub$dose, knots=knots[2:3], Bound=knots[c(1,4)])

# ALLOWING NON-LINEARITY IN BOTH FIXED AND RANDOM PARTS
nonlin <- mixmeta(logrr ~ 0 + spl, random= ~ 0 + spl|id, data=sub,
  method="ml", control=list(addSlist=addS))
summary(nonlin)

# SIMPLIFY THE MODEL BY ASSUMING LINEARITY IN THE RANDOM PART
nonlin2 <- update(nonlin, random= ~ 0 + dose|id)
summary(nonlin2)

# FIXED-EFFECTS MODEL (TRICK: random TO DEFINE THE GROUPING, THEN FIX IT TO 0)
nonlinfix <- mixmeta(logrr ~ 0 + spl, random= ~ 1|id, data=sub, method="ml",
  bscov="fixed", control=list(addSlist=addS, Psifix=0))
summary(nonlinfix)

# COMPARE THE MODELS
AIC(nocor, lin, nonlin, nonlin2, nonlinfix)

# PREDICT (RECREATE SPLINES FOR EASY CODING)
predlin <- exp(predict(lin, newdata=data.frame(dose=0:60), ci=TRUE))
newdata <- list(spl=ns(0:60, knots=knots[2:3], Bound=knots[c(1,4)]))
prednonlin <- exp(predict(nonlin, newdata=newdata, ci=TRUE))

# DISPLAY THE NON-LINEAR EFFECT
col1 <- do.call(rgb, c(as.list(col2rgb("blue") / 255), list(0.2)))
col2 <- do.call(rgb, c(as.list(col2rgb("green") / 255), list(0.2)))
plot(0:60, predlin[,1], type="l", ylim=c(0.85,1.9), ylab="RR",
  xlab="Alcohol intake (gr/day)", main="Dose-response")
polygon(c(0:60,60:0), c(predlin[,2], rev(predlin[,3])), col=col1, border=NA)
lines(0:60,prednonlin[,1], lty=5)
polygon(c(0:60,60:0), c(prednonlin[,2],rev(prednonlin[,3])), col=col2, border=NA)
}

\keyword{datasets}
