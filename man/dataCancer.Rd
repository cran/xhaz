% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataCancer.R
\docType{data}
\name{dataCancer}
\alias{dataCancer}
\title{Simulated data with cause death information with non comparability bias in term of individuals expected hazard}
\format{
This dataset contains the following variables:
\describe{
\item{obs_time}{Follow-up time (months)}
\item{obs_time_year}{Follow-up time (years)}
\item{event}{Vital status}
\item{age}{Age at diagnosis}
\item{agegrp}{"<30" , "30_60" and  ">=60" age groups }
\item{ageCentre}{centered age at diagnosis}
\item{sexx}{Sex(Female,Male).}
\item{immuno_trt}{Treatment group}
\item{year_date}{date of diagnosis.}
}
}
\usage{
data(dataCancer)
}
\description{
Simulated data
}
\examples{
data(dataCancer)
summary(dataCancer)
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
}
\keyword{datasets}
