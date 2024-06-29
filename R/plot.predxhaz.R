#' @title plots of excess hazard and net Survival from an
#'  \code{predxhaz} object
#'
#' @description Function to plot excess hazard or net survival
#'
#'
#' @param x An object of class predxhaz
#'
#' @param what allow to choose between excess hazard
#'      (\code{what="hazard"}) or net survival (\code{what="survival"}).
#'
#' @param ... additional arguments affecting the plot function
#'
#' @keywords plot.predxhaz
#'
#' @return The return of this function produce graphics of excess hazard or net
#'  survival, or time-dependent effects, when times.pts argument is provided
#'  in prediction call.
#'
#' @author Juste Goungounga, Robert Darlin Mba, Nathalie Graff\'eo and Roch Giorgi
#'
#' @importFrom graphics plot lines grid legend
#'
#' @export
#'
#'
#' @references Goungounga JA, Touraine C, Graff\'eo N, Giorgi R;
#' CENSUR working survival group. Correcting for misclassification
#' and selection effects in estimating net survival in clinical trials.
#'  BMC Med Res Methodol. 2019 May 16;19(1):104.
#'   doi: 10.1186/s12874-019-0747-3. PMID: 31096911; PMCID: PMC6524224.
#' (\href{https://pubmed.ncbi.nlm.nih.gov/31096911/}{PubMed})
#'
#' Touraine C, Graff\'eo N, Giorgi R; CENSUR working survival group.
#' More accurate cancer-related excess mortality through correcting
#' background mortality for extra variables.
#'  Stat Methods Med Res. 2020 Jan;29(1):122-136.
#'  doi: 10.1177/0962280218823234. Epub 2019 Jan 23. PMID: 30674229.
#'  (\href{https://pubmed.ncbi.nlm.nih.gov/30674229/}{PubMed})
#'
#'  Mba RD, Goungounga JA, Graff\'eo N, Giorgi R; CENSUR working survival group.
#'   Correcting inaccurate background mortality in excess hazard models
#'   through breakpoints. BMC Med Res Methodol. 2020 Oct 29;20(1):268.
#'   doi: 10.1186/s12874-020-01139-z. PMID: 33121436; PMCID: PMC7596976.
#'   (\href{https://pubmed.ncbi.nlm.nih.gov/33121436/}{PubMed})
#'
#' @examples
#' \donttest{
#' data("dataCancer")
#' # load the data set in the package
#' library("survival")
#' library("numDeriv")
#' library("survexp.fr")
#' data("simuData", package = "xhaz") # load the data sets 'simuData'
#'
#' #define the levels of variable sex
#'
#' # Esteve et al. model
#'
#'
#' fit.estv1 <- xhaz(formula = Surv(time_year, status) ~ agec + race,
#'                   data = simuData, ratetable = survexp.us,
#'                   interval = c(0, NA, NA, NA, NA, NA, max(simuData$time_year)),
#'                   rmap = list(age = 'age', sex = 'sex', year = 'date'),
#'                   baseline = "constant", pophaz = "classic")
#'
#'
#' predict_est <- predict(object = fit.estv1,
#'                        new.data = simuData,
#'                        times.pts = c(seq(0, 4, 0.1)),
#'                        baseline = TRUE)
#'
#' plot(predict_est, what = "survival",
#'      xlab = "time since diagnosis (year)",
#'      ylab = "net survival", ylim = c(0, 1))
#
#' data("dataCancer", package = "xhaz")   # load the data set in the package
#'
#' fit.phBS <- xhaz(
#'       formula = Surv(obs_time_year, event) ~ ageCentre + immuno_trt,
#'       data = dataCancer, ratetable = survexp.fr::survexp.fr,
#'       interval = c(0, NA, NA, max(dataCancer$obs_time_year)),
#'       rmap = list(age = 'age', sex = 'sexx', year = 'year_date'),
#'       baseline = "bsplines", pophaz  = "classic")
#'
#'
#' predict_mod1 <- predict(object = fit.phBS, new.data = dataCancer,
#'                         times.pts = c(seq(0, 10, 0.1)), baseline = FALSE)
#'
#' old.par <- par(no.readonly = TRUE)
#' par(mfrow = c(2, 1))
#'
#'
#' plot(predict_mod1, what = "survival",
#'      xlab = "time since diagnosis (year)",
#'      ylab = "net survival", ylim = c(0, 1))
#'
#' plot(predict_mod1, what = "hazard",
#'      xlab = "time since diagnosis (year)",
#'      ylab = "excess hazard")
#'
#' par(old.par)
#' }

plot.predxhaz <- function(x, what = "survival", ...){
   if (any(class(x) == "predxhaz")) {
    time <- sapply(1:length(x), function(i)unique(x[[i]]$times.pts))
    if (what == "survival") {
      survival <- sapply(1:length(x), function(i) mean(x[[i]]$survival))
      plot(time, survival, type = "l",...)
      grid()

    } else if (what == "hazard") {
      hazard <-
        sapply(1:length(x), function(i) {
          (sum(x[[i]]$hazard * x[[i]]$survival) / sum(x[[i]]$survival))
        })


      if (attr(x, "baseline") == "constant") {
        plot(time, hazard, type = "s",...)

        grid()

         }else {
        plot(time, hazard, type = "l",...)
        grid()
      }
    } else if (what == "beta") {
      stop("not yet implemented")
    }

  } else {
    stop("not yet implemented")
  }

}
