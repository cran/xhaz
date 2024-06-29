#' @title anova.mexhazLT function used for likelihood-ratio Test of two models
#' from mexhaz function
#'
#' @description  This function compute an analysis of deviance table for two
#' excess hazard models fitted using xhaz R package.
#'
#' @param object an object of class mexhazLT
#'
#' @param ... an object of class mexhazLT
#'
#' @param test a character string. The appropriate test is a likelihood-ratio
#' test, all other choices result in Not yet implemented test.
#'
#' @keywords anova.mexhazLT
#'
#' @note As expected, the comparison between two or more models by anova or more
#' excess hazard models will only be valid if they are fitted to the same
#' dataset, and if the compared models are nested. This may be a problem if
#' there are missing values.
#'
#'
#'
#' @return An object of class \code{anova} inheriting from class \code{matrix}.
#' The different columns contain respectively the degrees of freedom and the
#' log-likelihood values of the two nested models, the degree of freedom of the
#' chi-square statistic, the chi-square statistic and the p-value of the
#' likelihood ratio test.
#'
#'
#' @seealso \code{\link{xhaz}}, \code{\link{mexhazLT}}, \code{\link{AIC.mexhazLT}}
#'
#'
#' @author Juste Goungounga, Hadrien Charvat, Robert Darlin Mba, Nathalie Graff\'eo and Roch Giorgi
#'
#'
#' @references Goungounga JA, Touraine C, Graff\'eo N, Giorgi R;
#' CENSUR working survival group. Correcting for misclassification
#' and selection effects in estimating net survival in clinical trials.
#'  BMC Med Res Methodol. 2019 May 16;19(1):104.
#'   doi: 10.1186/s12874-019-0747-3. PMID: 31096911; PMCID: PMC6524224.
#' (\href{https://pubmed.ncbi.nlm.nih.gov/31096911/}{PubMed})
#'
#'  Goungounga, JA, Graff\'eo N, Charvat H, Giorgi R. “Correcting for
#'  heterogeneity and non-comparability bias in multicenter clinical trials
#'  with a rescaled random-effect excess hazard model.” Biometrical journal.
#'  Biometrische Zeitschrift vol. 65,4 (2023): e2100210.
#'  doi:10.1002/bimj.202100210.PMID: 36890623;
#'  (\href{https://pubmed.ncbi.nlm.nih.gov/36890623/}{PubMed})
#'
#'
#'
#' @examples
#' \donttest{
#' # load the data set in the package
#' library("survival")
#' library("numDeriv")
#' library("survexp.fr")
#'
#'
#' breast$sexe <- "female"
#'
#' fit.haz <- exphaz(
#'                   formula = Surv(temps, statut) ~ 1,
#'                   data = breast, ratetable = survexp.us,
#'                   only_ehazard = FALSE,
#'                   rmap = list(age = 'age', sex = 'sexe', year = 'date'))
#'
#' breast$expected <- fit.haz$ehazard
#' breast$expectedCum <- fit.haz$ehazardInt
#'
#' mod.bs3 <- mexhazLT(formula = Surv(temps, statut) ~ agecr + armt,
#'                   data = breast,
#'                   ratetable = survexp.us, degree = 3,
#'                   knots=quantile(breast[breast$statut==1,]$temps, probs=c(1:2/3)),
#'                   expected = "expected",expectedCum = "expectedCum",
#'                   base = "exp.bs", pophaz = "classic", random ="hosp")
#'
#' mod.bs3
#'
#' mod.bs4 <- mexhazLT(formula = Surv(temps, statut) ~ agecr + armt,
#'                   data = breast,
#'                   ratetable = survexp.us, degree = 3,
#'                   knots=quantile(breast[breast$statut==1,]$temps, probs=c(1:2/3)),
#'                   expected = "expected",expectedCum = "expectedCum",
#'                   base = "exp.bs", pophaz = "rescaled", random = "hosp")
#'
#' mod.bs4
#'
#'   anova(mod.bs3, mod.bs4)
#' }
#'
#' @importFrom stats printCoefmat
#' @export
anova.mexhazLT <- function(object, ...,
                           test = "LRT") {

if (test == "LRT") {

  if (!inherits(object, "mexhazLT"))
    stop("argument must be a xhaz fitted model")

  args <- list(...)

  if (length(args) >= 1 & any("mexhazLT" %in% unlist(lapply(1:length(args),
                                                         function(i) {
                                                           class(args[[i]])
                                                         })))) {
    nmodels <- length(unlist(lapply(1:length(args),
                                    function(i) {
                                      inherits(args[[i]], "mexhazLT")
                                    })))
  } else{
    nmodels <- 0
  }

  if (nmodels == 1) {
    object2 <- args[[1]]
  }else{
    stop("The anova function compare only two models")
  }

  if (!inherits(object2, "mexhazLT"))
    stop("argument must be a xhaz fitted model")

  if (length(object$loglik) > 1) {
    pvalue <- 1 - pchisq(2*(abs(object$loglik[2] - object2$loglik[2])),
                         df = abs(length(object$coefficients) - length(object2$coefficients)))
  } else if (length(object$loglik) == 1) {
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
    Chisq <- c(NA, round(abs(object2$loglik[2] - object$loglik[2]), 3))
  } else if (length(object$loglik) == 1) {
    loglik <- c(object$loglik[1], object2$loglik[1])
    dif.df <- c(NA, abs(length(object2$coef) - length(object$coef)))
    Chisq <- c(NA, round(abs(object2$loglik[1] - object$loglik[1]), 3))
  }


  p.value <- c(NA, round(pvalue, 10))
  x <- cbind(df, loglik, dif.df, Chisq, p.value)
  colnames(x) <- c("Model.df", "loglik",
                   "df", "Chisq", "Pr(>Chisq)")

  class(x) <- c("anova","matrix", "array" )


  printCoefmat(x,
               P.values = TRUE,
               digits = max(getOption("digits") - 2L, 3L),
               signif.stars = TRUE,
               na.print = "NA", has.Pvalue = TRUE)


}else {
  stop("Not yet implemented test!")
}
  invisible(x)
}



