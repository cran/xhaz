% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/exphaz.R
\name{exphaz}
\alias{exphaz}
\title{exphaz function}
\usage{
exphaz(
  formula = formula(data),
  data = sys.parent(),
  ratetable,
  rmap = list(age = NULL, sex = NULL, year = NULL),
  ratedata = sys.parent(),
  only_ehazard = TRUE,
  subset,
  na.action,
  scale = 365.2425
)
}
\arguments{
\item{formula}{a formula object of the \code{\link{Surv}} function with the
response on the left of a \code{~} operator and the terms on the right. The
response must be a survival object as returned by the \code{\link{Surv}}
function (\code{time} in first and \code{status} in second).}

\item{data}{a data frame in which to interpret the variables named in the
formula}

\item{ratetable}{a rate table stratified by \code{age}, \code{sex},
\code{year} (if missing, ratedata is used)}

\item{rmap}{a list that maps data set names to the ratetable names.}

\item{ratedata}{a data frame of the hazards mortality in general population.}

\item{only_ehazard}{a boolean argument (by default, \code{only_ehazard=TRUE}).
If \code{TRUE}, the cumulative population hazard is not provided.}

\item{subset}{an expression indicating which subset of the rows in data
should be used in the fit. All observations are included by default}

\item{na.action}{a missing-data filter function. The default is na.fail,
which returns an error if any missing values are found. An alternative is
na.exclude, which deletes observations that contain one or more missing
values.}

\item{scale}{a numeric argument specifying if the ratetable contains death
rates per day (default \code{scale = 365.2425}) or death rates per
year (\code{scale = 1}).}
}
\value{
An object of class \code{list} containing the following components:

\item{ehazard}{expected hazard calculated from the matching \code{ratetable}.}

\item{ehazardInt}{cumulative expected hazard calculated from the matching \code{ratetable}. if \code{only_ehazard=TRUE}, this quantity is not provided.}

\item{dateDiag}{date of diagnosis}
}
\description{
Calculate the expected hazard and survival.
}
\note{
\code{Time} is OBLIGATORY in YEARS.
}
\examples{

library(survexp.fr)
library(xhaz)
fit.haz <- exphaz(
                formula = Surv(obs_time_year, event) ~ 1,
                data = dataCancer,
                ratetable = survexp.fr, only_ehazard = TRUE,
                rmap = list(age = 'age', sex = 'sexx', year = 'year_date')
)

}
\references{
Goungounga JA, Touraine C, Graff??o N, Giorgi R;
CENSUR working survival group. Correcting for misclassification
and selection effects in estimating net survival in clinical trials.
BMC Med Res Methodol. 2019 May 16;19(1):104.
doi: 10.1186/s12874-019-0747-3. PMID: 31096911; PMCID: PMC6524224.
(\href{https://pubmed.ncbi.nlm.nih.gov/31096911/}{PubMed})

Touraine C, Graff??o N, Giorgi R; CENSUR working survival group.
More accurate cancer-related excess mortality through correcting
background mortality for extra variables.
Stat Methods Med Res. 2020 Jan;29(1):122-136.
doi: 10.1177/0962280218823234. Epub 2019 Jan 23. PMID: 30674229.
(\href{https://pubmed.ncbi.nlm.nih.gov/30674229/}{PubMed})

Mba RD, Goungounga JA, Graff??o N, Giorgi R; CENSUR working survival group.
Correcting inaccurate background mortality in excess hazard models
through breakpoints. BMC Med Res Methodol. 2020 Oct 29;20(1):268.
doi: 10.1186/s12874-020-01139-z. PMID: 33121436; PMCID: PMC7596976.
(\href{https://pubmed.ncbi.nlm.nih.gov/33121436/}{PubMed})
}
