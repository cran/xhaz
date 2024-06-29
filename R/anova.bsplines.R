#' @title anova.bsplines function used for likelihood-ratio Test of two models
#' from xhaz function
#'
#' @description  This function compute an analysis of deviance table for two
#' excess hazard models fitted using xhaz R package.
#'
#' @param object an object of class bsplines
#'
#' @param ... an object of class bsplines
#'
#' @param test a character string. The appropriate test is a likelihood-ratio
#' test, all other choices result in Not yet implemented test.
#'
#'
#' @keywords anova.bsplines
#'
#' @note As expected, the comparison between two or more models by anova or more
#' excess hazard models will only be valid if they are fitted to the same
#' dataset, and if the compared models are nested. This may be a problem if
#' there are missing values.
#'
#' @return An object of class \code{anova} inheriting from class \code{matrix}.
#' The different columns contain respectively the degrees of freedom and the
#' log-likelihood values of the two nested models, the degree of freedom of the
#' chi-square statistic, the chi-square statistic and the p-value of the
#' likelihood ratio test.
#'
#' @seealso \code{\link{xhaz}}, \code{\link{summary.bsplines}}, \code{\link{print.constant}}
#'
#'
#' @author Juste Goungounga, Robert Darlin Mba, Nathalie Graff\'eo and
#'  Roch Giorgi
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
#'
#'
#' @examples
#'
#' # load the data set in the package
#' \donttest{
#' library("survival")
#' library("numDeriv")
#' library("survexp.fr")
#' library("statmod")
#'
#' data("dataCancer", package = "xhaz")   # load the data set in the package
#'
#' fit.phBS <- xhaz(
#'       formula = Surv(obs_time_year, event) ~ ageCentre + immuno_trt,
#'       data = dataCancer,
#'       ratetable = survexp.fr::survexp.fr,
#'       interval = c(0, NA, NA, max(dataCancer$obs_time_year)),
#'       rmap = list(age = 'age', sex = 'sexx', year = 'year_date'),
#'       baseline = "bsplines", pophaz  = "classic")
#'
#'
#'
#' fit.nphBS <- xhaz(
#'       formula = Surv(obs_time_year, event) ~ ageCentre + qbs(immuno_trt),
#'       data = dataCancer,
#'       ratetable = survexp.fr::survexp.fr,
#'       interval = c(0, NA, NA, max(dataCancer$obs_time_year)),
#'       rmap = list(age = 'age', sex = 'sexx', year = 'year_date'),
#'       baseline = "bsplines", pophaz  = "classic")
#'
#' anova(fit.phBS, fit.nphBS)
#' }
#'
#' @importFrom stats printCoefmat
#' @export
anova.bsplines <- function(object, ..., test = "LRT") {
  if (test == "LRT") {
    if (!inherits(object, "bsplines"))
      stop("argument must be a xhaz fitted model")


    args <- list(...)

    if (length(args) >= 1 & any("bsplines" %in% unlist(lapply(1:length(args),
                                                              function(i) {
                                                                class(args[[i]])
                                                              })))) {
      nmodels <- length(unlist(lapply(1:length(args),
                                      function(i) {
                                        inherits(args[[i]], "bsplines")
                                      })))
    } else{
      nmodels <- 0
    }

    if (nmodels == 1) {
      object2 <- args[[1]]
    }else{
      stop("The anova function compare only two models")
    }


    if (!inherits(object2, c("bsplines", "constant")))
      stop("argument must be a xhaz fitted model")

    if (length(object$loglik) > 1) {
      pvalue <- 1 - pchisq(2*(abs(object$loglik[2] - object2$loglik[2])),
                           df = abs(length(object$coefficients) - length(object2$coefficients)))
    }else if (length(object$loglik) == 1) {
      pvalue <- 1 - pchisq(2*(abs(object$loglik[1] - object2$loglik[1])),
                           df = abs(length(object$coefficients) - length(object2$coefficients)))
    }

    cat("Assumption: Model 1 nested within Model 2
\n")
    cat("Likelihood ratio test\n")
    cat("Model 1: \n")
    print(object$call[2][[1]])
    cat("Model 2: \n")
    print(object2$call[2][[1]])
    df <- c(length(object$coef), length(object2$coef))
    if (length(object$loglik) > 1) {
      loglik <- c(object$loglik[2], object2$loglik[2])
      dif.df <- c(NA, abs(length(object2$coef) - length(object$coef)))
      Chisq <- c(NA, round(abs(object2$loglik[2] - object$loglik[2]),
                           3))
    }else if (length(object$loglik) == 1) {
      loglik <- c(object$loglik[1], object2$loglik[1])
      dif.df <- c(NA, abs(length(object2$coef) - length(object$coef)))
      Chisq <- c(NA, round(abs(object2$loglik[1] - object$loglik[1]),
                           3))
    }

    p.value <- c(NA, round(pvalue, 10))
    x <- cbind(df, loglik, dif.df, Chisq, p.value)
    colnames(x) <- c("Model.df", "loglik",
                     "Df", "Chisq", "Pr(>Chisq)")

    class(x) <- c("anova","matrix", "array" )


    printCoefmat(x,
                 P.values = TRUE,
                 digits = max(getOption("digits") - 2L, 3L),
                 signif.stars = TRUE,
                 na.print = "",
                 has.Pvalue = TRUE)

  }else {
    stop("Not yet implemented test!")
  }

  invisible(x)
}

