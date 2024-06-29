#' @title exphaz function
#'
#' @description Calculate the expected hazard and survival.
#'
#' @param formula a formula object of the \code{Surv} function with the
#' response on the left of a \code{~} operator and the terms on the right. The
#' response must be a survival object as returned by the \code{Surv}
#' function (\code{time} in first and \code{status} in second).
#' @note \code{Time} is OBLIGATORY in YEARS.
#'
#'
#' @param data a data frame in which to interpret the variables named in the
#' formula
#'
#' @param ratetable a rate table stratified by \code{age}, \code{sex},
#' \code{year} (if missing, ratedata is used)
#'
#' @param rmap a list that maps data set names to the ratetable names.
#'
#' @param ratedata a data frame of the hazards mortality in general population.
#'
#' @param only_ehazard a boolean argument (by default, \code{only_ehazard=TRUE}).
#' If \code{TRUE}, the cumulative population hazard is not provided.
#'
#' @param subset an expression indicating which subset of the rows in data
#' should be used in the fit. All observations are included by default
#'
#' @param na.action a missing data filter function. The default is na.fail,
#' which returns an error if any missing values are found. An alternative is
#' na.exclude, which deletes observations that contain one or more missing
#' values.
#'
#'
#' @param scale a numeric argument specifying by default \code{scale = 365.2425}
#' (or using the value corresponding to \code{attributes(ratetable)$cutpoints[[1]][2]}, often equal
#'  to 365.25) if the user wants to extract a yearly hazard rate, or \code{scale = 1} if he
#'  wants to extract a daily hazard rate from a ratetable containing daily hazard rates for
#'  a matched subject from the population, defined as \code{-log(1-q)/365.25}
#'  where \code{q} is the \code{1-}year probability of death.
#'
#' @return An object of class \code{list} containing the following components:
#'
#'
#' \item{ehazard}{expected hazard calculated from the matching \code{ratetable}.}
#'
#' \item{ehazardInt}{cumulative expected hazard calculated from the matching \code{ratetable}. if \code{only_ehazard=TRUE}, this quantity is not provided.}
#'
#' \item{dateDiag}{date of diagnosis}
#'
#' @references Goungounga JA, Touraine C, Graff\'eo N, Giorgi R;
#' CENSUR working survival group. Correcting for misclassification
#' and selection effects in estimating net survival in clinical trials.
#' BMC Med Res Methodol. 2019 May 16;19(1):104.
#' doi: 10.1186/s12874-019-0747-3. PMID: 31096911; PMCID: PMC6524224.
#' (\href{https://pubmed.ncbi.nlm.nih.gov/31096911/}{PubMed})
#'
#' Therneau, T. M., Grambsch, P. M., Therneau, T. M., & Grambsch, P. M. (2000).
#' Expected survival. Modeling survival data: extending the Cox model, 261-287.
#'
#' @examples
#'
#' library(survival)
#' library(survexp.fr)
#' library(xhaz)
#' fit.haz <- exphaz(
#'                 formula = Surv(obs_time_year, event) ~ 1,
#'                 data = dataCancer,
#'                 ratetable = survexp.fr, only_ehazard = TRUE,
#'                 rmap = list(age = 'age', sex = 'sexx', year = 'year_date')
#' )
#'
#' @export
exphaz <- function(formula = formula(data),
                 data = sys.parent(),
                 ratetable, rmap = list(age = NULL, sex = NULL, year = NULL),
                 ratedata = sys.parent(),
                 only_ehazard = TRUE,
                 subset,
                 na.action,
                 scale = 365.2425) {
  Call <- match.call()
  m <- match.call(expand.dots = FALSE)
  indx <- match(c("formula", "data", "subset", "na.action"),
                names(Call),
                nomatch = 0)
  if (indx[1] == 0)
    stop("A formula argument is required")
  temp <- Call[c(1, indx)]
  temp[[1]] <- as.name("model.frame")
  special <- c("strata")
  Terms <- if (missing(data)) {
    terms(formula, special)
  }
  else{
    terms(formula, special, data = data)
  }
  temp$formula <- Terms

  m <- eval(temp, sys.parent())
  if (missing(na.action))
    na.action <- NULL
  ehazardInt <- NULL
  # controls on data & ratetable parameters
  if (missing(ratedata) & missing(ratetable)) {
    stop("Missing rate table from general population.")
  }
  if (missing(data)) {
    stop("Missing data data frame in which to interpret
           the variables named in the formula.")
  } else{
    if (is.na(match(rmap$age, names(data))))
      stop("Must have informations for age on the data set.")
    if (is.na(match(rmap$sex, names(data))))
      stop("Must have informations for sex on the data set.")
    if (is.na(match(rmap$year, names(data))))
      stop("Must have informations for date on the data set.")
  }


  if (!missing(ratetable)) {
    if (is.ratetable(ratetable)) {
      varlist <- attr(ratetable, "dimid")
      if (is.null(varlist)) {
        varlist <- names(attr(ratetable, "dimnames"))
      }
      if (is.null(attributes(ratetable)$dimid)) {
        attributes(ratetable)$dimid <- varlist
      }
    }
    else{
      stop("Invalid rate table")
    }



    varsexID <- try(which(varlist == 'sex'))
    conditionVsex <- attr(ratetable, which = "dimnames")[[varsexID]]
    condsexdata <- unique(data[, rmap$sex])
    if (any(conditionVsex %in% condsexdata)) {
      conditionVsex <- condsexdata
    } else{
      stop(
        "Please check the matching between the levels of sex \nin the data.frame and in the ratetable used."
      )
    }



    if (!missing(rmap)) {
      rcall <- substitute(rmap)
      if (!is.call(rcall) || rcall[[1]] != as.name("list"))
        stop("Invalid rcall argument")
    }
    else
      rcall <- NULL

    temp01 <- match(names(rcall)[-1], varlist)
    if (any(is.na(temp01)))
      stop("Variable not found in the ratetable:",
           (names(rcall))[is.na(temp01)])

    temp02 <- match(as.vector(unlist(rmap)), names(data))
    if (any(is.na(temp02))) {
      stop("Variable not found in the data set:",
           (names(rcall))[is.na(temp02)])
    }

  }


  myvarnames <- colnames(model.matrix(Terms, m)[,-1, drop = FALSE])

  Y <- model.extract(m, "response")
  if (!inherits(Y, "Surv"))
    stop("Response must be a survival object.")

  strats <- attr(Terms, "specials")$strata
  dropx <- NULL

  attr(Terms, "intercept") <- 1

  if (length(dropx)) {
    X <- model.matrix(Terms[-dropx], m)[, -1, drop = FALSE]
  } else{
    X <- model.matrix(Terms, m)[, -1, drop = FALSE]
  }

  ###If there is a time-dependent covariate
  if (ncol(Y) == 2) {
    time <- Y[, 1]
    event <- Y[, 2]
  } else{
    time <- Y[, 2] - Y[, 1]
    event <- Y[, 3]
  }

  ageDiag <- data[, rmap$age]

if (missing(ratetable)) {
  exphaz <- exphaz_years(
    ageDiag = ageDiag,
    time = time,
    data = data,
    rmap = rmap,
    ratetable = ratetable,
    varlist = varlist,
    temp01 = temp01,
    scale = scale,
    pophaz = "rescaled",
    only_ehazard = only_ehazard
  )
  ehazard <- exphaz$ehazard
  ehazardInt <- try(exphaz$ehazardInt, TRUE)
} else{
  exphaz <- exphaz_years(
    ageDiag = ageDiag,
    time = time,
    data = data,
    rmap = rmap,
    ratetable = ratetable,
    varlist = varlist,
    temp01 = temp01,
    scale = scale,
    pophaz = "rescaled",
    only_ehazard = only_ehazard
  )
  if (only_ehazard) {
    ehazard <- exphaz$ehazard
    dateDiag <- exphaz$dateDiag

    return(list(ehazard = ehazard,
                dateDiag = dateDiag))
  }else {
    ehazard <- exphaz$ehazard
    ehazardInt <- exphaz$ehazardInt
    dateDiag <- exphaz$dateDiag

    return(list(ehazard = ehazard,
                ehazardInt = ehazardInt,
                dateDiag = dateDiag))
  }

}

}
