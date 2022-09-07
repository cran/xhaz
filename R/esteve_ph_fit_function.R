esteve.ph.fit <- function(x, y, ehazard, ehazardInt, int, covtest, bsplines,
                          init, control, event, Terms,strats, add.rmap,
                          add.rmap.cut, ageDiag, ageDC, optim, trace, speedy) {

  k <- length(int) - 1
  nvar <- ncol(x)
  nrowx <- nrow(x)

  nstrata <- ifelse(is.null(strats), 1, max(strats))
  attr(Terms, "nstrata") <- nstrata

  if (!is.null(add.rmap)) {
    if (!is.factor(add.rmap)) {
      stop("The alpha argument must be a factor.")
    } else{
      if (add.rmap.cut$breakpoint) {
        nalpha <- (length(add.rmap.cut$cut) + 1)*nlevels(add.rmap)

      }else{
        nalpha <- nlevels(add.rmap)
      }

    }
  } else{
    nalpha <- 0
  }


  #Initialization of theta0
  if (is.null(init)) {
    theta0 <- c(rep(0, nvar + k * nstrata + nalpha))
  } else{
    if (length(unlist(init)) != (nvar + k * nstrata + nalpha)) {
      stop("The number of initials values must the same as
             \nthe number of parameters to estimate.")
    } else{
      theta0 <- unlist(init)
      names(theta0) <- NULL
    }
  }

  #Initialization of tau
  tau <- c(rep(0, k * nstrata))
  names(tau) <- rep(sapply(1:k,
                           function(i, int1) {
                             paste("[", int1[i],
                                   "-",
                                   int1[i + 1],
                                   "[", sep = "")},
                           int1 = round(int, 2)), nstrata)

  if (nalpha) {
    if (nalpha == 1) {
      names(theta0) <- c(dimnames(x)[[2]], names(tau), "alpha")
    } else{
      if (add.rmap.cut$breakpoint) {


          break.levels <- levels(cut(ageDC, breaks = c(min(ageDC), add.rmap.cut$cut, max(ageDC))))



        names(theta0) <- c(dimnames(x)[[2]],
                           names(tau),
                           c(unlist(
                             lapply(1:nlevels(add.rmap),
                                    function(i)
                                      paste(
                                        paste0("alpha.",levels(add.rmap))[i],
                                        break.levels,
                                        # 1:(length(add.rmap.cut$cut) + 1),
                                        sep = "_")))))







      } else {
        names(theta0) <- c(dimnames(x)[[2]],
                           names(tau),
                           paste0('alpha.',
                                  levels(add.rmap)))
      }

    }
  } else{
    names(theta0) <- c(dimnames(x)[[2]], names(tau))
  }


  integr1 <- matrix(rep(0, nrowx * k), ncol = k)
  indic <- matrix(rep(0, nrowx * k), ncol = k)



  #The indicator function is different if a there is a time-dependent covariate
  if (ncol(y) == 2) {
    integr1 <- sapply(1:k, function(i, int, y)
      (int[i + 1] < y[, 1]) * (int[i + 1] - int[i]) +#/ 12 +
          (int[i] <= y[, 1] & y[, 1] <= int[i + 1]) * (y[, 1] - int[i]), #/ 12,
      int = int, y = y)
    indic <- sapply(1:k, function(i, int, y)
      (int[i] <= y[, 1] & y[, 1] < int[i + 1]),
      int = int, y = y)
  } else {
    integr1 <- sapply(1:k, function(i, int, y)
      ((y[, 1] <= int[i]) &
           (int[i + 1] < y[, 2])) * (int[i + 1] - int[i]) + #/ 12 +
          ((y[, 1] <= int[i]) & (int[i] <= y[, 2] &
                y[, 2] <= int[i + 1])) * (y[, 2] - int[i]) +#/ 12 +
          ((y[, 1] > int[i] & y[, 1] < int[i + 1]) &
             (int[i + 1] < y[, 2])) * (int[i + 1] - y[, 1]) +#/ 12 +
          ((y[, 1] > int[i] & y[, 1] < int[i + 1]) &
             (int[i] <= y[, 2] & y[, 2] <= int[i + 1])) *
          (y[, 2] - y[, 1]), #/ 12,
      int = int, y = y)

    indic <- sapply(1:k, function(i, int, y)
      (int[i] <= y[, 2] & y[, 2] < int[i + 1]),
      int = int,
      y = y)}

  #Full model
  if (optim) {
    Fmodel <- esteve.ph.optim.maxim(x = x, y = y,
                                    theta0 = theta0,
                                    nvar = nvar,
                                    k = k,
                                    indic = indic,
                                    event = event,
                                    integr1 = integr1,
                                    ehazard = ehazard,
                                    ehazardInt = ehazardInt,
                                    control.iter.max = control$iter.max,
                                    control.eps = control$eps,
                                    Terms = Terms,
                                    strats = strats,
                                    nstrata = nstrata,
                                    add.rmap = add.rmap,
                                    add.rmap.cut = add.rmap.cut,
                                    ageDiag = ageDiag,
                                    trace = trace,
                                    speedy = speedy)
    }
  else{
    Fmodel <- esteve.ph.maxim(x = x, y,
                              theta0 = theta0,
                              nvar = nvar,
                              k = k,
                              indic = indic,
                              event = event,
                              integr1 = integr1,
                              ehazard = ehazard,
                              ehazardInt = ehazardInt,
                              control.iter.max = control$iter.max,
                              control.eps = control$eps,
                              Terms = Terms,
                              strats = strats,
                              nstrata = nstrata,
                              add.rmap = add.rmap,
                              add.rmap.cut = add.rmap.cut,
                              ageDiag = ageDiag,
                              trace = trace)
    }
  #Null model: only coefficients for the baseline
  if (nvar > 0) {
    theta0N <- Fmodel$theta0[(nvar + 1):(nvar + k * nstrata + nalpha)]
    nvarN <- 0
    xN <- as.matrix(rep(1, nrow(x)), ncol = 1)
    if (optim) {
      Nmodel <- esteve.ph.optim.maxim(x = xN, y = y,
                                      theta0 = theta0N,
                                      nvar = nvarN,
                                      k = k,
                                      indic = indic,
                                      event = event,
                                      integr1 = integr1,
                                      ehazard = ehazard,
                                      ehazardInt = ehazardInt,
                                      control.iter.max = control$iter.max,
                                      control.eps = control$eps,
                                      Terms = Terms,
                                      strats = strats,
                                      nstrata = nstrata,
                                      add.rmap = add.rmap,
                                      add.rmap.cut = add.rmap.cut,
                                      ageDiag = ageDiag,
                                      trace = trace,
                                      speedy = speedy)

      } else{
      Nmodel <- esteve.ph.maxim(x = xN, y = y,
                                theta0 = theta0N,
                                nvar = nvarN,
                                k = k,
                                indic = indic,
                                event = event,
                                integr1 = integr1,
                                ehazard = ehazard,
                                ehazardInt = ehazardInt,
                                control.iter.max = control$iter.max,
                                control.eps = control$eps,
                                Terms = Terms,
                                strats = strats,
                                nstrata = nstrata,
                                add.rmap = add.rmap,
                                add.rmap.cut = add.rmap.cut,
                                ageDiag = ageDiag,
                                trace = trace)
    }
  } else{
    Nmodel <- Fmodel
  }

  if (sum(covtest) > 0) {
    cov.test <- covtest <- rep(TRUE, ncol(x))

    xT <- as.matrix(x[, grep("FALSE", as.character(covtest))])
    theta0T <- c(Fmodel$theta0[grep("FALSE", as.character(covtest))],
                 Fmodel$theta0[(nvar + 1):(nvar + k * nstrata)])
    nvarT <- ncol(xT)

    Tmodel <- esteve.ph.maxim(xT, y,
                              theta0T,
                              nvarT,
                              k,
                              indic,
                              event,
                              integr1,
                              ehazard,
                              ehazardInt,
                              control$iter.max,
                              control$eps,
                              Terms,
                              strats,
                              nstrata,
                              add.rmap,
                              add.rmap.cut,
                              ageDiag,
                              trace = trace
                              )


    dum <- c(rep(0, length(theta0)))
    names(dum) <- names(theta0)
    Tmodel$theta0 <- replace(dum,
                             grep("FALSE",
                                  as.character(is.na(
                                    match(names(theta0),
                                          names(Tmodel$theta0))))),
                             Tmodel$theta0)

    TF <- esteve.ph.maxim(x, y,
                          Tmodel$theta0,
                          nvar,
                          k,
                          indic,
                          event,
                          integr1,
                          ehazard,
                          ehazardInt,
                          control$iter.max,
                          control$eps,
                          Terms,
                          strats,
                          nstrata,
                          add.rmap,
                          add.rmap.cut,
                          ageDiag,
                          trace = trace)

    loglik.test <- -2 * (Tmodel$ll - Fmodel$ll)
    names(loglik.test) <- NULL
    wald.test <- t(Fmodel$theta0 - Tmodel$theta0) %*%
      (Fmodel$SD) %*% (Fmodel$theta0 - Tmodel$theta0)

    score.test <- t(TF$FD) %*% solve(TF$SD) %*% TF$FD
    Tmodel$wald.test <- wald.test
    Tmodel$score.test <- score.test
    Tmodel$loglik.test <- loglik.test
    var <- solve(Fmodel$SD)
    dimnames(var) <- list(names(theta0), names(theta0))
    list(
      coefficients = Fmodel$theta0,
      var = var,
      loglik = c(Nmodel$ll, Fmodel$ll),
      wald.test = Tmodel$wald.test,
      score.test = Tmodel$score.test,
      loglik.test = Tmodel$loglik.test,
      iterations = Fmodel$iter,
      cov.test = cov.test,
      cov.df = (ncol(x) - nvarT)
    )



  } else{
    cov.test <- FALSE
    var <- solve(Fmodel$SD)

    dimnames(var) <- list(names(theta0), names(theta0))
    if (optim) {

      list(
        coefficients = Fmodel$theta0,
        var = var,
        loglik = c(Nmodel$ll, Fmodel$ll),
        iterations = Fmodel$iter,
        cov.test = cov.test,
        message = Fmodel$message,
        convergence = Fmodel$convergence
      )
    } else{
      list(
        coefficients = Fmodel$theta0,
        var = var,
        loglik = c(Nmodel$ll, Fmodel$ll),
        iterations = Fmodel$iter,
        cov.test = cov.test
      )
    }
  }
  }

