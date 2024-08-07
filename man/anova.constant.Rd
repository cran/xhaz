% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/anova.constant.R
\name{anova.constant}
\alias{anova.constant}
\title{anova.constant function used for likelihood-ratio Test of two models
from xhaz function}
\usage{
\method{anova}{constant}(object, ..., test = "LRT")
}
\arguments{
\item{object}{an object of class constant}

\item{...}{an object of class constant}

\item{test}{a character string. The appropriate test is a likelihood-ratio
test, all other choices result in Not yet implemented test.}
}
\value{
An object of class \code{anova} inheriting from class \code{matrix}.
The different columns contain respectively the degrees of freedom and the
log-likelihood values of the two nested models, the degree of freedom of the
chi-square statistic, the chi-square statistic and the p-value of the
likelihood ratio test.
}
\description{
This function compute an analysis of deviance table for two
excess hazard models fitted using xhaz R package.
}
\note{
As expected, the comparison between two or more models by anova or more
excess hazard models will only be valid if they are fitted to the same
dataset, and if the compared models are nested. This may be a problem if
there are missing values.
}
\examples{
\donttest{
# load the data set in the package
library("survival")
library("numDeriv")
library("survexp.fr")


data("dataCancer")   # load the data set in the package

fit.ph <- xhaz(
      formula = Surv(obs_time_year, event) ~ ageCentre + immuno_trt,
      data = dataCancer,
      ratetable = survexp.fr::survexp.fr,
      interval = c(0, NA, NA, NA, max(dataCancer$obs_time_year)),
      rmap = list(age = 'age', sex = 'sexx', year = 'year_date'),
      baseline = "constant", pophaz  = "classic")



fit.ph2 <- xhaz(
      formula = Surv(obs_time_year, event) ~ ageCentre ,
      data = dataCancer,
      ratetable = survexp.fr::survexp.fr,
      interval = c(0, NA, NA, NA, max(dataCancer$obs_time_year)),
      rmap = list(age = 'age', sex = 'sexx', year = 'year_date'),
      baseline = "constant", pophaz  = "classic")

  anova(fit.ph2, fit.ph)
}

}
\references{
Goungounga JA, Touraine C, Graff\'eo N, Giorgi R;
CENSUR working survival group. Correcting for misclassification
and selection effects in estimating net survival in clinical trials.
BMC Med Res Methodol. 2019 May 16;19(1):104.
doi: 10.1186/s12874-019-0747-3. PMID: 31096911; PMCID: PMC6524224.
(\href{https://pubmed.ncbi.nlm.nih.gov/31096911/}{PubMed})

Touraine C, Graff\'eo N, Giorgi R; CENSUR working survival group.
More accurate cancer-related excess mortality through correcting
background mortality for extra variables.
Stat Methods Med Res. 2020 Jan;29(1):122-136.
doi: 10.1177/0962280218823234. Epub 2019 Jan 23. PMID: 30674229.
(\href{https://pubmed.ncbi.nlm.nih.gov/30674229/}{PubMed})

Mba RD, Goungounga JA, Graff\'eo N, Giorgi R; CENSUR working survival group.
Correcting inaccurate background mortality in excess hazard models
through breakpoints. BMC Med Res Methodol. 2020 Oct 29;20(1):268.
doi: 10.1186/s12874-020-01139-z. PMID: 33121436; PMCID: PMC7596976.
(\href{https://pubmed.ncbi.nlm.nih.gov/33121436/}{PubMed})

Giorgi R, Abrahamowicz M, Quantin C, Bolard P, Esteve J, Gouvernet J,
Faivre J. A relative survival regression model using B-spline functions
to model non-proportional hazards.
Statistics in Medicine 2003; 22: 2767-84.
(\href{https://pubmed.ncbi.nlm.nih.gov/12939785/}{PubMed})
}
\seealso{
\code{\link{xhaz}}, \code{\link{summary.bsplines}}, \code{\link{print.constant}}
}
\author{
Juste Goungounga, Robert Darlin Mba, Nathalie Graff\'eo and Roch Giorgi
}
\keyword{anova.constant}
