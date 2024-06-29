#' @import survival
#' @import stats
#' @import parallel
#' @import optimParallel

xhaz_split <- function(formula = formula(data),
                 data = sys.parent(),
                 ratetable, rmap = list(age = NULL, sex = NULL, year = NULL),
                 baseline = c("constant", "bsplines"),
                 pophaz = c("classic", "rescaled", "corrected"),
                 only_ehazard = FALSE,
                 add.rmap = NULL,
                 add.rmap.cut = list(breakpoint = FALSE,
                                     cut = c(70),
                                     probs = NULL,
                                     criterion = "BIC",
                                     print_stepwise = FALSE),
                 splitting = FALSE,
                 interval,
                 covtest,
                 ratedata = sys.parent(),
                 subset,
                 na.action,
                 init,
                 control = list(eps = 1e-4,
                                iter.max = 800,
                                level = 0.95),
                 optim = TRUE,
                 scale = 365.2425,
                 trace = 0,
                 speedy = FALSE,
                 nghq = 12, rcall, ...) {
  time_elapsed0 <- as.numeric(base::proc.time()[3])
  m <- match.call(expand.dots = FALSE)
  Call <- match.call()
  indx <- match(c("formula", "data", "subset", "na.action"),
                names(Call),
                nomatch = 0)


  if (indx[1] == 0)
    stop("A formula argument is required")
  temp <- Call[c(1, indx)]
  temp[[1L]] <- as.name("model.frame")


  special <- c("strata")
  Terms <- if (missing(data)) {
    terms(formula, special)
  }
  else{
    terms(formula, special, data = data)
  }
  temp$formula <- Terms
  m <- eval(temp, sys.parent())



  if (missing(na.action)) {
    na.action <- NULL
  } else if (length(attr(m, "na.action"))) {
    temp$na.action <- na.pass
    m <- eval(temp, sys.parent())
  }

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
    if (any(!conditionVsex %in% c('male', 'female'))) {
      conditionVsex <- c('male', 'female')[c(which(conditionVsex %in% c('male', 'female')))]
    }


    if (!missing(rmap)) {
      if (!splitting & missing(rcall)) {
        rcall <- substitute(rmap)
      } else if (!splitting & !missing(rcall)) {
        rmap <- eval(rmap)
      }
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


  if (pophaz == "corrected") {
    if (is.null(add.rmap.cut$breakpoint)) {
      stop("Missing breakpoint information")
    } else {
      if (add.rmap.cut$breakpoint == TRUE) {

        if (!is.na(add.rmap.cut$cut[1])) {

          if (min(add.rmap.cut$cut) < min(c(data[, rmap$age], data[, rmap$age] + data$time))) {
            if (max(add.rmap.cut$cut) <= max(c(data[, rmap$age], data[, rmap$age] + data$time))) {
              stop("Breakpoint(s) is (are) smaller than the minimum age")
            } else
              stop(
                "Breakpoint(s) is (are) smaller than the minimum age and breakpoint(s) greater than the maximum age"
              )
          } else{
            if (max(add.rmap.cut$cut) > max(c(data[, rmap$age], data[, rmap$age] + data$time)))
              stop("Breakpoint(s) is (are) greater than the maximum age")
          }

        }

      }
    }
  }



  if (control$iter.max < 0)
    stop("Invalid value for iterations.")
  if (control$eps <= 0)
    stop("Invalid convergence criteria.")
  if (control$level < 0 | control$level > 1)
    stop("Invalid value for the level of confidence interval.")

  if (missing(init))
    init <- NULL

  if (missing(interval))
    stop("Missing cutpoints definition for intervals.")
  if (!is.numeric(interval))
    stop("Wrong values for intervals. Must be numeric.")


  if (min(interval, na.rm = TRUE) != 0)
    stop("First interval must start at 0.")
  if (sum((interval < 0) * 1, na.rm = TRUE) > 0)
    stop("Negative value is not allowed for interval.")


  myvarnames <- colnames( model.matrix(Terms, m)[,-1, drop = FALSE])
  qbs_id <- which(stringr::str_detect(c(myvarnames),
                                      pattern = "qbs"))
  if (length(qbs_id) > 0) {

    if (length(interval) > 4)
      stop(
        "Interval must have 4 values using bsplines
          (2 internal knots plus '0' and the end of the study)."
      )
  }else{
    if (baseline == "bsplines") {
      if (length(interval) > 4)
        stop(
          "Interval must have 4 values using bsplines
          (2 internal knots plus '0' and the end of the study)."
        )
    }
  }

  Y <- model.extract(m, "response")
  if (!inherits(Y, "Surv"))
    stop("Response must be a survival object.")

  strats <- attr(Terms, "specials")$strata
  dropx <- NULL


  if (length(strats)) {
    if (length(qbs_id) > 0)
      stop("Strata function is not yet implemented for the B-splines model.")
    temp <- untangle.specials(Terms, "strata", 1)
    dropx <- c(dropx, temp$terms)
    if (length(temp$vars) == 1)
      strata.keep <- m[[temp$vars]]
    else
      strata.keep <- strata(m[, temp$vars], shortlabel = TRUE)
    strats <- as.numeric(strata.keep)
    attr(Terms, "nstrata") <- max(strats)
  }

  attr(Terms, "intercept") <- 1




  if (length(dropx)) {
    X <- model.matrix(Terms[-dropx], m)[, -1, drop = FALSE]
  } else{
    X <- model.matrix(Terms, m)[, -1, drop = FALSE]
  }


  if (length(qbs_id) > 0) {
    z_bsplines <- as.data.frame(
      model.matrix(Terms, m)[,-1, drop = FALSE][, c(qbs_id)])
    z_bsplines_names <- stringr::str_remove(myvarnames[c(qbs_id)],
                                            "qbs")
    colnames(z_bsplines) <-  gsub("\\(|\\)",
                                  "",
                                  as.character(z_bsplines_names))

    colnames(X)[c(qbs_id)] <- colnames(z_bsplines)
    z_bsplines <- as.matrix(z_bsplines)
    z_bsplines_vect <- rep(TRUE, ncol(z_bsplines))
    z_X_vect <- rep(FALSE, ncol(X))
    z_X_vect[c(qbs_id)] <- z_bsplines_vect
  }


  type <- attr(Y, "type")

  ###If there is a time-dependent covariate
  if (ncol(Y) == 2) {
    time <- Y[, 1]
    event <- Y[, 2]
  } else{
    time <- Y[, 2] - Y[, 1]
    event <- Y[, 3]
  }
  event[time > max(interval, na.rm = TRUE)] <- 0
  time[time > max(interval, na.rm = TRUE)] <- max(interval, na.rm = TRUE)
  if (length(qbs_id) > 0) {
    Y[, 1] <- time
  }
  ageDiag <- data[, rmap$age]
  ageDC <- ageDiag + time



  return(list(X = X,
              Y = Y,
              time = time,
              event = event,
              ageDC = ageDC,
              ageDiag = ageDiag))



}
