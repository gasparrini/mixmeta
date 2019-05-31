\name{dbs}
\alias{dbs}
\docType{data}

\title{ Deep-Brain Stimulation for Patients with Parkinson's Disease }

\description{
The dataset contains the data on 46 studies published between 1980 and 2004 that assessed the effect of deep-brain stimulation on the relief of symptoms of Parkinson's disease. The outcome is reported as a score motor function, defined with the Unified Parkinson's Disease Rating Scale (UPDRS-part III), with lower values indicating better prognosis. Changes in the score were measured at 3, 6, 12 months and long-term after the implantation of the stimulator.
}

\usage{dbs}

\format{
  A data frame with 68 observations on the following 12 variables:
  \itemize{
    \item{\code{author}: label identifying the study.}
    \item{\code{year}: year of publication.}
    \item{\code{eff_month3,  var_month3}: point estimate and variance of the change in the score at 3 months.}
    \item{\code{eff_month6,  var_month6}: point estimate and variance of the change in the score at 6 months.}
    \item{\code{eff_month12,  var_month12}: point estimate and variance of the change in the score at 12 months.}
    \item{\code{eff_long,  var_long}: point estimate and variance of the change in the score in the long term.}
    \item{\code{duration}: average disease duration (years).}
    \item{\code{baseline}: average baseline score of the patients.}
  }
}

\details{
The data are stored in a \emph{wide} format, with each record belonging to a single study and different variables providing estimates of the outcome at different times. Each study report results at one or multiple times, with the remaining times set to missing. See the dataset \link{gliomas} for an example of similar dataset stored in \emph{long} format.
}

\source{
Ishak KJ, et al (2007). Meta-analysis of effect sizes reported at multiple time points using general linear mixed model. \emph{Clinical Trials}. \bold{4}(5):525--39.

Sera F, Gasparrini A (2019). An extended mixed-effects framework for meta-analysis.\emph{Statistics in Medicine}. In submission.
}

\note{
The data provide an example of application of longitudinal meta-analysis, with repeated measurements of the effect size taken at various time point within each study. This requires a modelling structure that accounts for both within and between-study correlations of repeated measurements. In this case, the analysis is performed in the wide-format dataset using a multivariate meta-analysis. However, a long format is better suited for longitudinal meta-analysis, as it is applicable even when estimates are reported at different times in each study (see the examples in the help page of the dataset \code{\link{gliomas}}). Results can be compared with those reported Ishak and colleagues (2007).
}

\examples{
### REPRODUCE THE RESULTS IN ISHAK ET AL (2007), TABLES 1 AND 2

# CREATE THE OUTCOME AND WITHIN-STUDY MATRICES (THE LATTER WITHOUT CORRELATION)
y <- as.matrix(dbs[1:4*2+1])
S <- as.matrix(dbs[1:4*2+2])

# INDEPENDENT RANDOM EFFECTS (TABLE 1, FIRST MODEL)
mv1 <- mixmeta(y ~ 1, S, bscov="diag", data=dbs)
print(summary(mv1), digits=1, report="var")

# HETEROGENEOUS AR1 RANDOM-EFFECTS (TABLE 1, THIRD MODEL)
mv3 <- mixmeta(y ~ 1, S, bscov="har1", data=dbs)
print(summary(mv3), digits=1, report="var")

# BUILD THE LIST HETEROGENEOUS AR1 WITHIN-STUDY ERRORS (CORRELATION AT 0.97)
cormat <- 0.97^abs(col(matrix(1,4,4)) - row(col(matrix(1,4,4))))
addS <- lapply(seq(nrow(S)), function(i) inputcov(sqrt(S[i,]), cormat))
addS <- lapply(addS, function(x) x[apply(!is.na(x),1,any), apply(!is.na(x),2,any)])

# ADD HAR1 WITHIN-STUDY ERRORS (TABLE 1, FOURTH MODEL) USING addSlist
mv4 <- mixmeta(y ~ 1, bscov="har1", data=dbs, control=list(addSlist=addS))
print(summary(mv4), digits=1, report="var")


### USE A LONG FORMAT, AS MORE FLEXIBLE AND ALLOWS MORE COMPLEX MODELS

# RESHAPE THE DATASET
long <- reshape(dbs, direction="long", idvar="author", v.names=c("eff","var"),
  varying=list(1:4*2+1, 1:4*2+2))

# RE-RUN THE LAST (FOURTH) MODEL
mv4b <- mixmeta(eff ~ factor(time) - 1, random = ~ factor(time) -1 | author,
  bscov="har1", data=long, control=list(addSlist=addS))
print(summary(mv4b), digits=1, report="var")

# COMMON RANDOM EFFECTS (TABLE 1, SECOND MODEL)
mv2 <- mixmeta(eff ~ factor(time) - 1, var, random = ~ factor(time) -1 | author,
  bscov="id", data=long)
print(summary(mv2), digits=1, report="var")

# FOURTH MODEL WITH ADDITIONAL CENTERED META-PREDICTORS (TABLE 2)
mv4plus <- mixmeta(eff ~ factor(time) - 1 + I(duration-14) + I(baseline-52),
  random = ~ factor(time) -1 | author, bscov="har1", data=long,
  control=list(addSlist=addS))
print(summary(mv4plus), digits=1, report="var")


### SEE help(gliomas) FOR A COMPLEMENTARY EXAMPLE
}

\keyword{datasets}
