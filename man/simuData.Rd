% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simuData.R
\docType{data}
\name{simuData}
\alias{simuData}
\title{Simulated data with cause death information in long term follow-up setting without non comparability bias in term of individuals expected hazard}
\format{
This dataset contains the following variables:
\describe{
\item{age}{Age at diagnosis}
\item{agec}{Centered age}
\item{sex}{Sex(Female,Male)}
\item{race}{Race}
\item{date}{date of diagnosis.}
\item{time}{Follow-up time (months)}
\item{time_year}{Follow-up time (years)}
\item{status}{Vital status}

}
}
\usage{
data(simuData)
}
\description{
Simulated data
}
\examples{
data(simuData)
summary(simuData)
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
}
\keyword{datasets}
