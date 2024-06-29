#' colorectum cancer data with multiple events
#'
#' multiple events data
#'
#'
#' @docType data
#'
#' @usage data(dataCancer)
#'
#' @format This dataset contains the following variables:
#' \describe{
#'  \item{id}{patient IDs.}
#'  \item{sex}{gender with 1 for male and 2 for female.}
#'  \item{sexe}{gender male and female.}
#'  \item{age}{Age at diagnosis}
#'  \item{stage}{lower to higher stage 1, 2, 3}
#'  \item{time}{time-to-events (local or distant recurrence or death)}
#'  \item{status}{ 0 : no event; 1: local or distant recurrence or death}
#'  \item{event}{1: local recurrence; 2: distant recurrence; 3:death}
#'  \item{date_diag}{date of diagnosis.}
#' }
#'
#'
#' @keywords datasets
#'
#' @references Touraine C, Graff\'eo N, Giorgi R; CENSUR working survival group.
#' More accurate cancer-related excess mortality through correcting
#' background mortality for extra variables.
#' Stat Methods Med Res. 2020 Jan;29(1):122-136.
#' doi: 10.1177/0962280218823234. Epub 2019 Jan 23. PMID: 30674229.
#' (\href{https://pubmed.ncbi.nlm.nih.gov/30674229/}{PubMed})
#'
#'
#' @examples
#' data(ccr.mevents)
#' summary(ccr.mevents)
"ccr.mevents"
