\name{gliomas}
\alias{gliomas}
\docType{data}

\title{ Randomized Trials on Therapies for Malignant Gliomas }

\description{
The dataset contains the data on 17 randomized controlled trials comparing post-operative radiation therapy plus chemotherapy versus radiation therapy alone in patients with malignant gliomas. The outcome of interest is the probability of survival along time, measured as odds ratio at 6, 12, 18, and 24 months.
}

\usage{gliomas}

\format{
  A data frame with 68 observations on the following 8 variables:
  \itemize{
    \item{\code{study}: number identifying the trial.}
    \item{\code{time}: time (months) since the start of the treatment at which survival status is assessed.}
    \item{\code{ntreat, streat}: number of total patients at the beginning of the study and surviving patients at specific times, respectively, in the treatment group (radiation therapy plus chemotherapy).}
    \item{\code{dcontr, ncontr}: number of total patients at the beginning of the study and surviving patients at specific times, respectively, in the control group (radiation alone).}
    \item{\code{logOR, varOR}: log-odds ratio of survival between treatment and control groups.}
  }
}

\details{
The data are stored in a \emph{long} format, with each record providing the estimate at a single time and each study providing multiple records. There were missing data for study 17 at months 6 and 18. There were no survivors in the control group at month 24 for studies 3 and 10, although this still allows computation of the OR. See the dataset \code{\link{dbs}} for an example of similar dataset stored in \emph{wide} format.

The log-odds ratio is computed empirically as \eqn{\log(s_t \times (n_t-s_t) / ((n_c-s_c) \times s_c)}{log(st x (nt-st) / ((nc-sc) x sc)}. Its variance is simply computed as \eqn{1/s_t + 1/(n_t-s_t) + 1/(n_c-s_c) + 1/s_c}{1/st + 1/(nt-st) + 1/(nc-sc) + 1/sc}.
}

\source{
Musekiwa A, et al (2012). Meta-analysis of effect sizes reported at multiple time points using general linear mixed model. \emph{Plos One}. \bold{11}(10):e0164898.

Trikalinos TA, Olkin I (2012). Meta-analysis of effect sizes reported at multiple time points: a multivariate approach. \emph{Clinical Trials}. \bold{9}(5):610--620.

Sera F, Armstrong B, Blangiardo M, Gasparrini A (2019). An extended mixed-effects framework for meta-analysis.\emph{Statistics in Medicine}. 2019;38(29):5429-5444. [Freely available \href{http://www.ag-myresearch.com/2019_sera_statmed.html}{\bold{here}}].
}

\note{
The data provide an example of application of longitudinal meta-analysis, with repeated measurements of the effect size taken at various time points within each study. This requires a modelling structure that accounts for both within and between-study correlations of repeated measurements. In this case, the same analysis can be performed in a wide-format dataset using a multivariate meta-analysis (see the examples in the help page of the dataset \code{\link{dbs}}). However, the long format is better suited for longitudinal meta-analysis, as it is applicable even when estimates are reported at different times in each study. Results can be compared with those reported Musekiwa and colleagues (2016). The same dataset was also used by Trikalinos and Olkin (2012), using a similar modelling scheme.
}

\examples{
### REPRODUCE THE RESULTS IN MUSEKIWA ET AL (2012), TABLES 3 AND 4

# INDEPENDENT RANDOM EFFECTS, NO WITHIN-STUDY CORRELATION (MODEL 1)
mod1 <- mixmeta(logOR~0+factor(time), S=logORvar, random=~0+factor(time)|study,
  bscov="diag", data=gliomas)
print(summary(mod1), digits=2, report="var")

# COMPOUND-SYMMETRY RANDOM EFFECTS, NO WITHIN-STUDY CORRELATION (MODEL 2)
# NB: THIS REQUIRES A TWO-LEVEL MODEL WITH THE INNER-LEVEL VARIANCE FIXED TO 0
unit <- factor(seq(nrow(gliomas)))
mod2 <- mixmeta(logOR~0+factor(time), S=logORvar, random=~1|study/unit,
  bscov=c("unstr","fixed"), data=gliomas, control=list(Psifix=list(unit=0)))
print(summary(mod2), digits=2, report="var")

# BUILD THE HETEROGENEOUS AR1 WITHIN-STUDY ERRORS (CORRELATION AT 0.61)
cormat <- 0.61^abs(col(matrix(1,4,4)) - row(col(matrix(1,4,4))))
addS <- lapply(split(sqrt(gliomas$logORvar), gliomas$study), inputcov, cormat)
addS <- lapply(addS, function(x) x[apply(!is.na(x),1,any), apply(!is.na(x),2,any)])

# INDEPENDENT RANDOM EFFECTS, HAR1 WITHIN-STUDY CORRELATION (MODEL 4)
mod4 <- mixmeta(logOR~0+factor(time), random=~0+factor(time)|study,
  bscov="diag", data=gliomas, control=list(addSlist=addS))
print(summary(mod4), digits=2, report="var")

# UNSTRUCTURED RANDOM EFFECTS, HAR1 WITHIN-STUDY CORRELATION (MODEL 6)
mod6 <- update(mod4, bscov="unstr")
print(summary(mod6), digits=2, report="var")

# COMPARE THE FIT
AIC(mod1, mod2, mod4, mod6)

### MORE FLEXIBLE MODELLING OF RANDOM EFFECTS

# RE-RUN BEST FITTING MODEL WITH ML (ALLOWS TESTING OF FIXED EFFECTS)
mod4ml <- update(mod4, method="ml")

# RANDOM-SLOPE MODEL WITH TIME AS CONTINUOUS (CENTERED IN random)
modnew <- mixmeta(logOR~time, random=~I(time-15)|study, bscov="diag",
  method="ml", data=gliomas, control=list(addSlist=addS, maxiter=200))
print(summary(modnew), digits=2, report="var")

# COMPARE
AIC(mod4ml, modnew)


### SEE help(dbs) FOR A COMPLEMENTARY EXAMPLE
}

\keyword{datasets}
