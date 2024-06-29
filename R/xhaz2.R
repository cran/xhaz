xhaz2 <- function(formula = formula,
                  data = data,
                  ratetable = ratetable,
                  rmap = rmap,
                  baseline = baseline,
                  pophaz = pophaz,
                  only_ehazard = only_ehazard,
                  add.rmap = add.rmap,
                  add.rmap.cut = add.rmap.cut,
                  splitting  = splitting,
                  interval = interval,
                  ratedata = ratedata,
                  subset = subset,
                  na.action = na.action,
                  init = init,
                  control = control,
                  optim = optim,
                  scale = scale,
                  trace = trace,
                  speedy = speedy,
                  nghq = nghq,
                  m_int = m_int,
                  rcall = rcall,
                  method = method,
                  ...) {
  time_elapsed0 <- as.numeric(base::proc.time()[3])


  Call <- match.call()

  m <- match.call(expand.dots = FALSE)

  indx <- match(c("formula", "data"),
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


  if (missing(subset)) {
    subset <- NULL
  }
  if (missing(na.action)) {
    na.action <- "na.omit"
    m <- eval(temp, sys.parent())

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
    if (is.na(match(rmap$age, names(data)))) {
      stop("Must have informations for age on the data set.")
    }
    if (is.na(match(rmap$sex, names(data)))) {
      stop("Must have informations for sex on the data set.")
    }
    if (is.na(match(rmap$year, names(data)))) {
      stop("Must have informations for date on the data set.")
    }
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
      condition2 <-
        add.rmap.cut$breakpoint == TRUE &
        is.na(add.rmap.cut$cut[1]) & !is.null(add.rmap.cut$probs)

      if ((!splitting & missing(rcall)) & (condition2)) {
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
           rmap[[which(is.na(temp02))]])
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


  myvarnames <- colnames(model.matrix(Terms, m)[,-1, drop = FALSE])
  qbs_id <- which(stringr::str_detect(c(myvarnames),
                                      pattern = "qbs"))
  if (length(qbs_id) > 0) {
    if (length(interval) > 4)
      stop(
        "Interval must have 4 values using bsplines
          (2 internal knots plus '0' and the end of the study)."
      )
  } else{
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
    z_bsplines <-
      as.data.frame(model.matrix(Terms, m)[, -1, drop = FALSE][, c(qbs_id)])
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
    covtest <- z_X_vect
  } else{
    covtest <- rep(FALSE, ncol(X))
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
  time[time > max(interval, na.rm = TRUE)] <-
    max(interval, na.rm = TRUE)
  if (length(qbs_id) > 0) {
    Y[, 1] <- time
  }

  if (is.null(data$break_interval)) {
    ageDiag <- data[, rmap$age]
    ageDC <- ageDiag + time
  } else{
    ageDiag <- data$tstart
    ageDC <- data$tstop
  }


  pophaz <- match.arg(pophaz, c("classic", "rescaled", "corrected"))

  if (pophaz == "corrected") {
    if (!is.null(add.rmap)) {
      add.rmap.var <- add.rmap
      add.rmap <- data[, add.rmap]
    } else{
      stop("Additional demographic variable must be specified")
    }
  } else{
    if (pophaz == "rescaled") {
      if (!is.null(add.rmap)) {
        stop("Additional demographic variable is not required")
      } else{
        add.rmap <- as.factor(rep(1, nrow(data)))
      }
    }

    if (pophaz == "classic") {
      if (!is.null(add.rmap)) {
        stop("Additional demographic variable is not required")
      }
    }
  }

  if (only_ehazard == TRUE & pophaz != "classic") {
    stop("cumulative expected hazard if also required for this type of model")
  }


  #
  condition0 <- add.rmap.cut$breakpoint == FALSE
  condition1 <-
    add.rmap.cut$breakpoint == TRUE &
    !is.na(add.rmap.cut$cut[1]) & is.null(add.rmap.cut$probs)
  condition2 <-
    add.rmap.cut$breakpoint == TRUE &
    is.na(add.rmap.cut$cut[1]) & !is.null(add.rmap.cut$probs)
  if (!is.null(data$break_interval)) {
    if (missing(ratetable)) {
      exphaz <- exphaz_years(
        ageDiag = data$tstart,
        time = time,
        data = data,
        rmap = rmap,
        ratetable = ratetable,
        varlist = varlist,
        temp01 = temp01,
        scale = scale,
        pophaz = pophaz,
        add.rmap = add.rmap,
        only_ehazard = only_ehazard
      )
      ehazard <- exphaz$ehazard
      ehazardInt <- try(exphaz$ehazardInt, TRUE)
    } else{
      exphaz <- exphaz_years(
        ageDiag = data$tstart,
        time = time,
        data = data,
        rmap = rmap,
        ratetable = ratetable,
        varlist = varlist,
        temp01 = temp01,
        scale = scale,
        pophaz = pophaz,
        only_ehazard = only_ehazard
      )
      ehazard <- exphaz$ehazard
      ehazardInt <- exphaz$ehazardInt
      dateDiag <- exphaz$dateDiag
    }
  } else {
    if (missing(ratetable)) {
      exphaz <- exphaz_years(
        ageDiag = ageDiag,
        time = time,
        data = data,
        rmap = rmap,
        ratetable = ratetable,
        ratedata = ratedata,
        varlist = varlist,
        temp01 = temp01,
        scale = scale,
        pophaz = pophaz,
        add.rmap = add.rmap,
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
        pophaz = pophaz,
        only_ehazard = only_ehazard
      )
      ehazard <- exphaz$ehazard
      ehazardInt <- exphaz$ehazardInt
      dateDiag <- exphaz$dateDiag
    }

  }



  if (sum(is.na(interval)) > 0) {
    n.cut <- sum(is.na(interval))
    q.values <- cumsum(rep(1 / (n.cut + 1), n.cut))

    if (baseline == "bsplines" & (n.cut != 2)) {
      if (n.cut != 3) {
        q.values <- c(0.05, 0.95)
      }
      else {
        stop("Must have 2 internal knots using bsplines.")
      }
    }


    l.cut <- quantile(time[which(event %in% 1)], q.values)
    names(l.cut) <- NULL
    interval <- c(min(interval, na.rm = TRUE),
                  l.cut,
                  max(interval, na.rm = TRUE))
  }

  if ((length(interval) - 1) != sum(sapply(1:(length(interval) - 1),
                                           function(i, interval)
                                             (interval[i + 1] > interval[i]),
                                           interval = interval)))
    stop("Interval values are not in ascending order.")


  if (!missing(covtest)) {
    if ((sum(covtest) == ncol(X)) && (length(qbs_id) == 0))
      stop(
        "Do not use 'covtest' for this hypothesis.
          \nLikelihood ratio test of the full versus null model
          \nis always provided."
      )

    if (length(covtest) != ncol(X))
      stop(
        "Number of arguments of 'covtest' must be the same
          \nas the number of fitted binaries covariates or
          \nas the number of levels if data type is factor."
      )
  } else{
    covtest <- c(rep(FALSE, ncol(X)))
  }



  if (length(qbs_id) > 0) {
    if ((length(z_X_vect) != ncol(X)) ||
        (!is.logical(z_X_vect)) || (sum(is.na(z_X_vect)) > 0))

      stop(
        "Invalid values for 'qbs()':
          \nmust be well specified for covariable(s) used in the formula."
      )

    if (ncol(Y) > 2)
      stop(
        "Time-dependent covariate not yet implemented for
          \nnon-proportional hazards situation."
      )

    for (i in 1:length(z_X_vect))
      if ((z_X_vect[i] == FALSE) && (covtest[i] == TRUE) == TRUE)
        stop("You mustn't test a PH effect (covtest=TRUE) for
               \na PH covariate (z_X_vect=FALSE)!")
  } else{
    z_X_vect <- covtest <- rep(FALSE, ncol(X))
    if ((length(z_X_vect) != ncol(X)) ||
        (!is.logical(z_X_vect)) || (sum(is.na(z_X_vect)) > 0))

      stop("You mustn't test a PH effect (covtest=TRUE) for
               \na PH covariate (z_X_vect=FALSE)!")

    if (ncol(Y) > 2)
      stop(
        "Time-dependent covariate not yet implemented for
          \nnon-proportional hazards situation."
      )

    for (i in 1:length(z_X_vect))
      if ((z_X_vect[i] == FALSE) && (covtest[i] == TRUE) == TRUE)
        stop("You mustn't test a PH effect (covtest=TRUE) for
               \na PH covariate (z_X_vect=FALSE)!")

  }


  baseline <-  match.arg(baseline, c("constant", "bsplines"))

  if (baseline == "constant") {
    if (add.rmap.cut$breakpoint == FALSE) {
      fitter <- get("esteve.ph.fit")

      fit <- fitter(
        X,
        Y,
        ehazard,
        ehazardInt,
        int = interval,
        covtest,
        bsplines =  z_X_vect,
        init,
        control,
        event,
        Terms,
        strats,
        add.rmap,
        add.rmap.cut,
        ageDiag,
        ageDC,
        optim,
        trace,
        speedy,
        method
      )
    } else if (add.rmap.cut$breakpoint == TRUE &
               !is.na(add.rmap.cut$cut[1]) &
               is.null(add.rmap.cut$probs)) {
      fitter <- get("esteve.ph.fit")

      fit <- fitter(
        X,
        Y,
        ehazard,
        ehazardInt,
        int = interval,
        covtest,
        bsplines =  z_X_vect,
        init,
        control,
        event,
        Terms,
        strats,
        add.rmap,
        add.rmap.cut,
        ageDiag,
        ageDC,
        optim,
        trace,
        speedy,
        method
      )
    } else if (add.rmap.cut$breakpoint == TRUE &
               is.na(add.rmap.cut$cut[1]) &
               !is.null(add.rmap.cut$probs)) {
      fitter <- get("esteve.ph.fit")
      if (splitting) {
        nbreak <- length(add.rmap.cut$cut)
        allpos_break <-
          with(data, quantile(ageDC[event == 1], probs = c(add.rmap.cut$probs)))
        cuted <- gtools::permutations(n = length(allpos_break),
                                      r = nbreak,
                                      v = allpos_break)

        if (nbreak > 1) {
          cut2 <- unique(t(sapply(1:nrow(cuted), function(i)
            sort(cuted[i,]))))
        } else{
          cut2 <-
            unique(matrix(sapply(1:nrow(cuted), function(i)
              sort(cuted[i,])),
              ncol = 1))
        }


        nmodels <- nrow(cut2)


        tofit <- lapply(1:nmodels, function(i) {
          add.rmap.cut$cut <- cut2[i, ]
          newdata2 <- tosplit(
            formula = formula,

            add.rmap.cut = add.rmap.cut,
            data = data,
            rmap = rmap,
            interval = interval,
            subset
          )
          data <- newdata2$tdata2

          if (is.null(data$break_interval)) {
            ageDiag <- data[, rmap$age]
            ageDC <- ageDiag + time
          } else if (!is.null(data$break_interval)) {
            ageDiag <- data$tstart
            ageDC <- data$tstop
            time <- with(data, c(tstop - tstart))


            add.rmap <- data[, add.rmap.var]

          }


          if (!survival::is.ratetable(ratetable)) {
            exphaz2 <- exphaz_years(
              ageDiag = ageDiag,
              time = time,
              data = data,
              rmap = rmap,
              ratetable = ratetable,
              varlist = varlist,
              temp01 = temp01,
              scale = scale,
              pophaz = pophaz,
              add.rmap = add.rmap,
              only_ehazard = only_ehazard
            )
            ehazard2 <- exphaz2$ehazard
            ehazardInt2 <- try(exphaz2$ehazardInt, TRUE)
          } else{
            exphaz2 <- exphaz_years(
              ageDiag = ageDiag,
              time = time,
              data = data,
              rmap = rmap,
              ratetable = ratetable,
              varlist = varlist,
              temp01 = temp01,
              scale = scale,
              pophaz = pophaz,
              only_ehazard = only_ehazard
            )
            ehazard2 <- data$ehazard2 <- exphaz2$ehazard
            ehazardInt2 <- data$ehazardInt2 <- exphaz2$ehazardInt
            dateDiag2 <- data$dateDiag2 <- exphaz2$dateDiag
          }

          newfit <- xhaz_split(
            formula = formula,
            data = data,
            ratetable = ratetable,
            rmap = rmap,
            baseline  = baseline,
            pophaz = pophaz,
            only_ehazard = only_ehazard,
            add.rmap = add.rmap,
            add.rmap.cut = add.rmap.cut,
            splitting  = splitting,
            interval = interval,
            covtest = covtest,
            init = init,
            control = control,
            optim = optim,
            scale = scale ,
            trace = trace,
            speedy = speedy,
            nghq,
            rcall = rcall,
            ...
          )

          X <- newfit$X
          Y <- newfit$Y
          event <- newfit$event
          ageDC <- newfit$ageDC
          ageDiag <- newfit$ageDiag
          testM(
            X,
            Y,
            ehazard = ehazard2,
            ehazardInt = ehazardInt2,
            int = interval,
            covtest,
            bsplines =  z_X_vect,
            init,
            control,
            event,
            Terms,
            strats,
            add.rmap,
            add.rmap.cut,
            ageDiag = ageDiag,
            ageDC = ageDC,
            optim,
            trace,
            speedy,
            data
          )

        })

        if (length(which(stringr::str_detect(
          names(unlist(add.rmap.cut)), "print_stepwise"
        ))) > 0) {
          if (add.rmap.cut$print_stepwise) {
            sapply(1:length(tofit),
                   function(i) {
                     cat("Model:", i, "\n")
                     if (length(which(stringr::str_detect(
                       names(tofit[[i]]), "coefficients"
                     ))) > 0)  {
                       tofit[[i]]$n <- nrow(Y)

                       tofit[[i]]$level <- control$level
                       tofit[[i]]$interval <- interval
                       tofit[[i]]$n.events <- sum(event, na.rm = TRUE)
                       tofit[[i]]$formula <- as.vector(attr(Terms, "formula"))
                       tofit[[i]]$call <- m_int
                       tofit[[i]]$varcov <- tofit[[i]]$var
                       tofit[[i]][["var"]] <- NULL
                       tofit[[i]]$pophaz <- pophaz
                       tofit[[i]]$baseline <- baseline
                       tofit[[i]]$add.rmap <- add.rmap
                       tofit[[i]]$add.rmap.cut  <- add.rmap.cut
                       if (!splitting) {
                         tofit[[i]]$terms <- Terms
                         tofit[[i]]$assign <- attr(X, "assign")
                       }
                       oldClass(tofit[[i]]) <- "constant"
                     }

                     if (length(which(stringr::str_detect(
                       names(tofit[[i]]), "coefficients"
                     ))) > 0)  {
                       # xhaz:::print.constant(tofit[[i]])

                       print.constant(tofit[[i]])
                     } else {
                       cat("Model", i, ": No convergence \n")
                     }
                     cat("\n")
                   })
            cat("\n")
          }
        }



        allAIC <- suppressWarnings(sapply(1:length(tofit), function(i)
            as.numeric(try(tofit[[i]]$AIC, TRUE))))
        allBIC <- suppressWarnings(sapply(1:length(tofit), function(i)
            as.numeric(try(tofit[[i]]$BIC, TRUE))))
        browser()
        if (isTRUE(which(names(add.rmap.cut) %in% "criterion") > 0)) {
          if (add.rmap.cut$criterion == "AIC") {
            if (any(!is.na(allAIC))) {
              fit <- tofit[[which.min(allAIC)]]
              fit$add.rmap.cut$cut <- c(cut2[which.min(allAIC),])
               }else {
              stop("\nNo convergence for any model")
            }
          } else if (add.rmap.cut$criterion == "BIC") {
           if (any(!is.na(allBIC))) {
             fit <- tofit[[which.min(allBIC)]]
             fit$add.rmap.cut$cut <- c(cut2[which.min(allBIC),])
           }else {
             stop("\nNo convergence for any model")
           }

          }

        } else{
          fit <- tofit[[which.min(allBIC)]]
          fit$add.rmap.cut$cut <- c(cut2[which.min(allBIC),])
        }


        fit$data <- data

      } else{
        nbreak <- length(add.rmap.cut$cut)

        age_time <- ageDiag + time
        allpos_break <-
          with(data, quantile(age_time[event == 1], probs = c(add.rmap.cut$probs)))
        cuted <- gtools::permutations(n = length(allpos_break),
                                      r = nbreak,
                                      v = allpos_break)

        if (nbreak > 1) {
          cut2 <- unique(t(sapply(1:nrow(cuted), function(i)
            sort(cuted[i,]))))
        } else{
          cut2 <-
            unique(matrix(sapply(1:nrow(cuted), function(i)
              sort(cuted[i,])),
              ncol = 1))
        }


        nmodels <- nrow(cut2)

        tofit <- lapply(1:nmodels, function(i) {
          add.rmap.cut$cut <- cut2[i, ]
          testM(
            X,
            Y,
            ehazard,
            ehazardInt,
            int = interval,
            covtest,
            bsplines =  z_X_vect,
            init,
            control,
            event,
            Terms,
            strats,
            add.rmap,
            add.rmap.cut,
            ageDiag,
            ageDC,
            optim,
            trace,
            speedy,
            data
          )

        })

        if (length(which(stringr::str_detect(
          names(unlist(add.rmap.cut)), "print_stepwise"
        ))) > 0) {
          if (add.rmap.cut$print_stepwise) {
            sapply(1:length(tofit),
                   function(i) {
                     cat("Model:", i, "\n")
                     print(tofit[[i]])
                     cat("\n")
                   })
            cat("\n")
          }
        }

        allAIC <-
          suppressWarnings(sapply(1:length(tofit), function(i)
            as.numeric(try(tofit[[i]]$AIC, TRUE))))
        allBIC <-
          suppressWarnings(sapply(1:length(tofit), function(i)
            as.numeric(try(tofit[[i]]$BIC, TRUE))))
        if (which.min(allAIC) < 1) {
          stop("no convergence with the proposed breakpoints")
        }

        if (add.rmap.cut$criterion == "AIC") {
          fit <- tofit[[which.min(allAIC)]]
          fit$add.rmap.cut$cut <- c(cut2[which.min(allAIC),])
        } else if (add.rmap.cut$criterion == "BIC") {
          fit <- tofit[[which.min(allBIC)]]
          fit$add.rmap.cut$cut <- c(cut2[which.min(allBIC),])
        }
      }
    }


    oldClass(fit) <- "constant"
  }
  else {
    fitter <- get("giorgi.tdph.fit")
    fit <- fitter(
      X,
      Y,
      ehazard,
      ehazardInt,
      int = interval,
      covtest,
      bsplines = z_X_vect,
      init,
      control,
      event,
      Terms,
      strats,
      add.rmap,
      add.rmap.cut,
      ageDiag,
      ageDC,
      optim,
      trace,
      speedy,
      nghq,
      method
    )
    oldClass(fit) <- "bsplines"
    fit$z_bsplines <- z_X_vect
  }
  time_elapsed1 <- as.numeric(base::proc.time()[3])

  if (add.rmap.cut$breakpoint == TRUE &
      !is.na(add.rmap.cut$cut[1])) {
    fit$break.levels <-
      levels(cut(ageDC, breaks = c(
        min(ageDC), add.rmap.cut$cut, max(ageDC)
      )))
  } else if (add.rmap.cut$breakpoint == TRUE &
             is.na(add.rmap.cut$cut[1])) {
    fit$break.levels <-
      levels(cut(ageDC, breaks = c(
        min(ageDC), fit$add.rmap.cut$cut, max(ageDC)
      )))

  }




  fit$level <- control$level
  fit$interval <- interval
  fit$na.action <- na.action
  fit$n <- nrow(Y)
  fit$n.events <- sum(event, na.rm = TRUE)
  fit$formula <- as.vector(attr(Terms, "formula"))
  fit$call <- m_int
  fit$varcov <- fit$var
  fit[["var"]] <- NULL
  fit$pophaz <- pophaz
  fit$baseline <- baseline
  fit$add.rmap <- add.rmap
  fit$ehazard <- ehazard
  fit$ehazardInt <- ehazardInt
  fit$add.rmap.cut  <- add.rmap.cut
  fit$time_elapsed <- time_elapsed1 - time_elapsed0


  if (!splitting) {
    fit$data <- data
    fit$terms <- Terms
    fit$assign <- attr(X, "assign")

  }
  return(fit)
}
