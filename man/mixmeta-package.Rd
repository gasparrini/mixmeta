\name{mixmeta-package}
\alias{mixmeta-package}
\docType{package}

\title{
An Extended Mixed-Effects Framework for Meta-Analysis
}

\description{
The package \pkg{mixmeta} consists of a collection of functions to perform various meta-analytical models in \R through a unified mixed-effects framework, including standard univariate fixed and random-effects meta-analysis and meta-regression, and non-standard extensions such as multivariate, multilevel, longitudinal, and dose-response models.
}

\section{Modelling framework}{
Standard applications of meta-analysis amount to the pooling of estimates of a single effect size, here defined generally as outcome, collected as unique observations in a set of independent studies, together with a measure of uncertainty (usually standard errors). Fixed-effects models do not assume heterogeneity across studies, and the estimates are conditional on the set of studies collected in the meta-analysis. Random-effects meta-analysis, instead, allows a degree of heterogeneity among studies, assuming the (true but unobserved) study-specific outcomes as randomly sampled from a (usually hypothetical) population of studies. Meta-regression extends both fixed and random-effects methods by allowing the pooled outcome to depend on study-level meta-predictors.

However, this traditional setting can be limited for many modern applications of meta-analysis. For instance, studies can provide estimates of different outcomes. Alternatively, studies can collect multiple estimates of the same outcome, either longitudinally or referring to different groups or levels of a continuous variable. Similarly, studies can be clustered, or being characterized by a hierarchical structure (i.e., by country). In all these instances, the key assumption of independence across estimates is not met, and basic models must be extended to consider potentially complex correlation structures within and between studies. This leads to extension to multivariate, multilevel, longitudinal, or dose-response models for meta-analysis, among others.

A unified modelling framework can be defined by casting the meta-analytical problem as a linear mixed model. In general terms, we assume that there is a set of \eqn{n} observations of \eqn{k} different outcomes, representing \emph{units} of analysis aggregated in \eqn{i = 1,\dots,m} \emph{groups} that are considered independent. An extended random-effects meta-regression model for the \eqn{\boldsymbol{y}_i}{yi} outcomes in group \eqn{i} can be generally written as:

\deqn{\boldsymbol{y}_i = \boldsymbol{X}_i\boldsymbol{\beta} +  \boldsymbol{Z}_i\boldsymbol{b} + \boldsymbol{\epsilon}_i}{yi = Xi\beta + Zib + \epsiloni}
\deqn{\boldsymbol{b} \sim \mbox{N}(\boldsymbol{0},\boldsymbol{\Psi}) , \boldsymbol{\epsilon}_i \sim \mbox{N}(\boldsymbol{0},\boldsymbol{S}_i)}{b ~ N(0,\Psi) , \epsiloni ~ N(0,Si)}

Here, \eqn{\boldsymbol{X}_i\boldsymbol{\beta}}{Xi\beta} defines the fixed effects that represent the population-averaged outcomes in terms of \eqn{p} unit-level predictors in the design matrix \eqn{\boldsymbol{X}_i}{Xi} with fixed-effects coefficients \eqn{\boldsymbol{\beta}}{\beta}. The random part of the model \eqn{\boldsymbol{Z}_i\boldsymbol{b}}{Zib} describes the deviation from the population averages in terms of \eqn{q} unit-level predictors in the design matrix \eqn{\boldsymbol{Z}_i}{Zi} and random-effects coefficients \eqn{\boldsymbol{b}}{b}. The marginal (co)variance matrix \eqn{\boldsymbol{\Sigma}_i = \boldsymbol{Z}\boldsymbol{\Psi}\boldsymbol{Z}^t + \boldsymbol{S}_i}{\Sigmai = Z\PsiZt + Si} is given by the sum of within (assumed known) and between-group contributions, defined by (co)variance matrices \eqn{\boldsymbol{S}_i}{Si} and \eqn{\boldsymbol{\Psi}}{\Psi}, respectively.

All the models mentioned above, and other extensions, can be described as special cases of this unified framework. Specifically, in the standard random-effects univariate meta-analysis or meta-regression, each group represents a study with a single observation (\eqn{n=m}), where \eqn{\boldsymbol{Z_i=1}}{Zi=1} (\eqn{q=1}), and \eqn{\boldsymbol{y}_i}{yi}, \eqn{\boldsymbol{S}_i}{Si} and \eqn{\boldsymbol{\Psi}}{\Psi} are scalars. In fixed-effects models, \eqn{\boldsymbol{\Psi}}{\Psi} and \eqn{\boldsymbol{Z_i}}{Zi} do not exist. In multivariate models, the \eqn{k}-dimensional \eqn{\boldsymbol{y}_i}{yi} represents the different outcomes from study \eqn{i}, \eqn{\boldsymbol{X}_i}{Xi} is Kronecker-expanded to \eqn{k \times kp}{k x kp}, and \eqn{\boldsymbol{S}_i}{Si} and \eqn{\boldsymbol{\Psi}}{\Psi} are \eqn{k \times k}{k x k} matrices representing within and between-study correlations among outcomes, respectively. In multilevel models, where additional inner levels of grouping exist within each of the \eqn{m} outer-level groups, \eqn{q} is the sum of level-specific meta-predictors, while \eqn{\boldsymbol{\Psi}}{\Psi} and \eqn{\boldsymbol{Z}_i}{Zi} have a block-diagonal and column-binded (and expanded) forms, respectively, with each part referring to a different level. In longitudinal and dose-response models, repeated measures are accomodated in a similar way through random-effects grouping.
}

\section{Estimation methods}{
The aim is to estimate the \eqn{kp} coefficients \eqn{\boldsymbol{\beta}}{\beta} and, for random-effects models, the components of the between-group (co)variance matrix \eqn{\boldsymbol{\Psi}}{\Psi}. The parameters for the random part depend on the number of random-effects levels, and for each of them, on the number of random-effects meta-predictors and the structure of the related part of the (co)variance matrix, with a maximum of \eqn{kq(kq+1)/2} for single-level unstructured \eqn{\boldsymbol{\Psi}}{\Psi}.

Different estimators are implemented in the package \pkg{mixmeta}. The options available in the current version are:

\itemize{
  \item \bold{Fixed-effects}
  \item \bold{Maximum likelihood} (ML)
  \item \bold{Restricted maximum likelihood} (REML)
  \item \bold{Method of moments}
  \item \bold{Variance components}
}

The fixed-effects model is fitted through generalized least squares (GLS), assuming the (co)variance structure, composed by the within-study errors only, as completely known. Likelihood-based random-effects estimators, ML and REML, represent the most comprehensive implementation of the modelling framework, and allow the specification all the various models described in the previous section through a flexible definition of the random-effects structure. They rely on two alternative iterative optimization procedures, based on Newton-type and (restricted) iterative generalized least squares (IGLS and RIGLS) algorithms, respectively. Estimators based on semiparametric alternatives such as the non-iterative method of moments or the iterative variance components are also included, although they are only available for models with a basic random-effects structure. Further details on estimation methods are given in the related help pages.
}

\section{Functions included in the package}{
The main function in the package is \code{\link{mixmeta}}, which performs the various models illustrated above. This function resembles standard regression functions in \R, and specifies the model through regression formulae for fixed and random-effects (see \code{\link{mixmetaFormula}}). The function returns a list object of class \code{"mixmeta"} (see \code{\link{mixmetaObject}}).

The estimation is carried out internally through \code{\link{mixmeta.fit}}, a wrapper which accepts data in a specific format, then prepares the various data components and calls ad hoc estimation functions for fitting the models. Specifically, \code{\link{mixmeta.fixed}} is applied for fixed-effects models, while estimators for random-effects models are implemented in the functions \code{\link{mixmeta.ml}} and \code{\link{mixmeta.reml}} for (restricted) maximum likelihood, \code{\link{mixmeta.mm}} for the method of moments, and \code{\link{mixmeta.vc}} for variance components. For likelihood-based methods, alternative iterative optimizations methods are provided in two sets of functions implementing \code{\link[=ml.newton]{Newton-type}} and \code{\link[=ml.igls]{(R)IGLS}} algorithms used for maximizing the (restricted) likelihood. The former method applies specific \code{\link[=ml.loglik.fn]{likelihood functions}}. Various types of likelihood-based models are defined by separate regression formulae for fixed and random-effects (see \code{\link{mixmetaFormula}}). Specific \code{\link[=mixmetaCovStruct]{(co)variance structures}} for the between-group random effects at single or multiple levels are available. Fitting parameter options are set by \code{\link{mixmeta.control}}.

Method functions are available for objects of class \code{"mixmeta"} (see \code{\link{mixmetaObject}} for a complete list). The method \code{\link[=summary.mixmeta]{summary}} produces a list of class \code{"summary.mixmeta"} for summarizing the fit of the model and providing additional results. The method function \code{\link[=predict.mixmeta]{predict}} computes predicted values, optionally for a set of new values of the predictors. \code{\link[=blup.mixmeta]{blup}} gives the (empirical) best linear unbiased predictions for the set of studies used for estimation. Other default or specific method functions for regression can be used on objects of class \code{"mixmeta"}, such as \code{\link{fitted}} and \code{\link{residuals}}, \code{\link[=logLik.mixmeta]{logLik}}, \code{\link{AIC}} and \code{\link{BIC}}, or \code{\link{drop1}} and \code{\link{add1}}, among others.

Methods for \code{\link[=model.frame.mixmeta]{model.frame}}, \code{\link[=model.matrix.mixmeta]{model.matrix}}, and \code{\link[=terms.mixmeta]{terms}} are used to extract or construct the model frame, the design matrix, or the terms of the regression meta-analytical model, respectively. These specific methods for objects of class \code{"mixmeta"} are needed to appropriately deal with missing values and to account for model frames that include terms for both the fixed and random parts. In particular, methods for \code{\link[=na.omit.data.frame.mixmeta]{na.omit}} and \code{\link[=na.exclude.data.frame.mixmeta]{na.exclude}} are used to handle correctly missing values.

Simulations can be produced using the function \code{\link{mixmetaSim}} and the method function \code{\link[=simulate.mixmeta]{simulate}}, which return one or multiple sets of simulated outcomes for a group of studies. The function \code{\link{inputna}} and \code{\link{inputcov}} are used internally to augment the missing data values and to input missing correlations, respectively.

The method function \code{\link{qtest.mixmeta}} (producing an object with class of the same name) performs the (multivariate) Cochran Q test for (residual) heterogeneity. For multivariate models, the function returns both an overall estimate and those for each single outcome. The generic method function is \code{\link{qtest}}.

Printing functions for the objects of classes defined above are also provided. Other functions are used internally in the source code, and not exported in the namespace. For users interested in getting into details of the package structure, these functions can be displayed using the triple colon ('\code{:::}') operator. For instance, \code{mixmeta:::glsfit} displays the code of the function \code{glsfit}. Also, some comments are added in the original source code.
}

\section{Datasets and applications}{

The package includes several datasets used for applications of the extended meta-analytical framework. The related help pages provide examples of specific models, and fully demonstrate the flexibility of the extended meta-analytical framework. In particular:

\itemize{
  \item \bold{Standard meta-analysis} is illustrated using the dataset \code{\link{bcg}}, including examples of meta-regression.
  \item \bold{Multivariate meta-analysis} is performed using various datasets, including bivariate models (\code{\link{berkey98}}, \code{\link{hyp}}, \code{\link{p53}}) and multivariate models with three or more outcomes (\code{\link{fibrinogen}} and \code{\link{hsls}}). The examples describe also how to deal with missing data or missing within-group correlations, and multivariate meta-regression.
  \item \bold{Network meta-analysis} is shown in the dataset \code{\link{smoking}}. The examples illustrate an indirect mixed-treatment comparison including consistency and inconsistency models.
  \item \bold{Multilevel meta-analysis} is displayed in the examples of the datasets \code{\link{school}} and \code{\link{thrombolytic}}, and include data with multiple nested levels of grouping and/or repeated measures within each group.
  \item \bold{Dose-response meta-analysis} is illustrated in the dataset \code{\link{alcohol}}, using the recently proposed one-stage approach.
  \item \bold{Longitudinal meta-analysis} is performed using the datasets \code{\link{dbs}} and \code{\link{gliomas}}. The two sets of examples present different cases using data in wide and long format, respectively.
}
}

\section{Additional information}{
The \pkg{mixmeta} package is available on the Comprehensive R Archive Network (CRAN), with info at the related web page (\href{https://CRAN.R-project.org/package=mixmeta}{CRAN.R-project.org/package=mixmeta}). A development website is available on GitHub (\href{https://github.com/gasparrini/mixmeta}{github.com/gasparrini/mixmeta}). General information on the development and applications of this extended meta-analytical modelling framework, together with an updated version of the \R{} scripts for running the examples in published papers, can be found in GitHub (\href{https://github.com/gasparrini}{github.com/gasparrini}) or at the personal web page of the package maintainer (\href{http://www.ag-myresearch.com}{www.ag-myresearch.com}).

The package \pkg{mixmeta} is an extension of the package \pkg{mvmeta}, previously developed to perform multivariate meta-analytical models. The latter now depends on the former, and while both are still maintained, users are encouraged to switch to \pkg{mixmeta} as it represents a more general and updated option. A list of changes included in the current and previous versions of \pkg{mixmeta} can be found by typing:

\code{news(package="mixmeta")}

Use \code{citation("mixmeta")} to cite this package.
}

\note{
The package \pkg{mixmeta} provides a unified modelling framework to perform standard and non-standard meta-analytical models. However, some of these can also be fitted using routines available in other \R packages.

For instance, many packages such as \pkg{metafor}, \pkg{meta}, \pkg{rmeta} provide a more exhaustive and efficient set of methods for standard univariate meta-analysis and meta-regression, including a wide range or functions for specific plots and statistical tests.

Specific modelling extensions are also provided by other packages. For example, multivariate or multilevel models can be also be fitted using functions in \pkg{metafor} and \pkg{metaSEM}, while dose-response meta-analysis and meta-analysis of diagnostic measures can be performed using \pkg{dosresmeta} and \pkg{mada}, respectively.

See the CRAN Task View \href{https://cran.r-project.org/view=MetaAnalysis}{Meta-Analysis} for a comprehensive illustration of methods available in various \R packages.
}

\section{Warnings}{
This release of the package \pkg{mixmeta} has been tested with different simulated and real datasets. The functions generally perform well under several scenarios, and comparisons with alternative software implementations show good agreement. However, bugs and bad performance under un-tested conditions may not be excluded. Please report any error or unexpected behaviour to the e-mail address below.
}

\author{
Antonio Gasparrini and Francesco Sera

Maintainer: Antonio Gasparrini <\email{antonio.gasparrini@lshtm.ac.uk}>
}

\references{
Sera F, Gasparrini A (2019). An extended mixed-effects framework for meta-analysis.\emph{Statistics in Medicine}. DOI: 10.1002/sim.8362 [Freely available \href{http://www.ag-myresearch.com/2019_sera_statmed.html}{\bold{here}}].

Gasparrini A, Armstrong B, Kenward MG (2012). Multivariate meta-analysis for non-linear and other multi-parameter associations. \emph{Statistics in Medicine}. \bold{31}(29):3821--3839. [Freely available \href{http://www.ag-myresearch.com/2012_gasparrini_statmed.html}{\bold{here}}].

Pinheiro JC and Bates DM (2000). \emph{Mixed-Effects Models in S and S-PLUS}. New York, Springer Verlag.

Lindstrom MJ and Bates DM (1988). Newton-Raphson and EM algorithms for linear mixed-effects models for repeated-measures data. \emph{Journal of the American Statistical Association}. \bold{83}(404):1014--1022.

Goldstein H (1986). Multilevel mixed linear model analysis using iterative generalized least squares. \emph{Biometrika}. \bold{73}(1):43--56.

Goldstein H (1992). Efficient computational procedures for the estimation of parameters in multilevel models based on iterative generalized least squares. \emph{Computational Statistics \& Data Analysis}. \bold{13}(1):63--71.

Stram DO (1996). Meta-analysis of published data using a linear mixed-effects model. \emph{Biometrics}. \bold{52}(2):536--544.

Stevens JR, Taylor AM. Hierarchical dependence in meta-analysis. \emph{Journal of Educational and Behavioral Statistics}. \bold{34}(1):46--73.

Jackson D, Riley R, White IR (2011). Multivariate meta-analysis: Potential and promise. \emph{Statistics in Medicine}. \bold{30}(20);2481--2498.

Goldstein H, et al (2000). Meta-analysis using multilevel models with an application to the study of class size effects. \emph{Journal of the Royal Statistical Society: Series C (Applied Statistics)}. \bold{49}(3):399--412.

Crippa A, et al (2019). One-stage dose-response meta-analysis for aggregated data. \emph{Statistical Methods in Medical Research}.  \bold{28}(5):1579--1596.

Ishak KJ, Platt RW, Joseph L, et al (2007). Meta-analysis of longitudinal studies. \emph{Clinical Trials}. \bold{4}(5):525--539.
}

\keyword{package}
