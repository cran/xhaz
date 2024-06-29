#' Simulated data with cause death information with non comparability bias in term of individuals expected hazard
#'
#' Simulated data
#'
#'
#' @docType data
#'
#' @usage data(rescaledData)
#'
#' @format This dataset contains the following variables:
#' \describe{
#'  \item{time}{Follow-up time (months)}
#'  \item{status}{Vital status}
#'  \item{age}{Age at diagnosis}
#'  \item{age.c}{Centred age}
#'  \item{sex}{Sex(Female,Male)}
#'  \item{hormTh}{Treatment group variable}
#'  \item{date}{date of diagnosis}
#'
#' }
#'
#'
#' @keywords datasets
#'
#' @references Goungounga JA, Touraine C, Graff\'eo N, Giorgi R;
#' CENSUR working survival group. Correcting for misclassification
#' and selection effects in estimating net survival in clinical trials.
#' BMC Med Res Methodol. 2019 May 16;19(1):104.
#' doi: 10.1186/s12874-019-0747-3. PMID: 31096911; PMCID: PMC6524224.
#' (\href{https://pubmed.ncbi.nlm.nih.gov/31096911/}{PubMed})
#'
#'
#' @examples
#' data(rescaledData)
#' summary(rescaledData)
"rescaledData"
