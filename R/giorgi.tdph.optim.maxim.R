#' @import statmod
#' @import numDeriv
giorgi.tdph.optim.maxim <- function(x,
                                    theta0,
                                    nTD,
                                    nPH,
                                    event,
                                    ehazard,
                                    ehazardInt,
                                    P,
                                    p,
                                    k,
                                    nrowx,
                                    IntGL,
                                    cpti,
                                    cptj,
                                    cst,
                                    boundmin,
                                    boundmax,
                                    Int.chazard,
                                    Int.FDbase,
                                    Int.SDbase,
                                    int,
                                    control.iter.max,
                                    control.eps,
                                    Terms,
                                    add.rmap,
                                    add.rmap.cut,
                                    trace,
                                    speedy = FALSE,
                                    nghq = nghq,
                                    pophaz = pophaz, method = method) {
  if (is.null(trace)) {
    trace <- 0
  }
  alpha <- NULL
  no <- length(add.rmap)
  k <- 3
  nrowx <- nrow(x)
  ehazard <- c(ehazard)
  event <- c(event)
  cst <- 1
  cpti <- 1
  cptj <- 1
  rrBetaPH <- 1
  rrBetaTD <- 1
  int.chazard <- matrix(0, nrowx, k)
  int.FDbase <- matrix(0, nrowx, 5 * k)
  int.SDbase <- matrix(0, 5 * nrowx, 5 * k)

  #Loop for integral calculation
  for (cst in 1:k) {
    int.chazard[, cst] <-
      IntGL(Int.chazard, t(matrix(c(
        boundmax[, cst], boundmin[, cst]
      ), ncol = 2)), cst, cpti, cptj, theta0[1:(5 + 5 * nTD)],
      x, nTD, p, nghq = nghq)
  }

  name.bs <- dimnames(x)[[2]]
  name <- c(rep(0, 5 + 5 * (nTD) + nPH))
  name[1:5] <- c("qbs base ( 1 )",
                 "qbs base ( 2 )",
                 "qbs base ( 3 )",
                 "qbs base ( 4 )",
                 "qbs base ( 5 )")

  if (!is.null(add.rmap)) {
    nalpha <- nlevels(add.rmap)
    rescale <- TRUE
  } else{
    if (pophaz == "classic") {
      rescale <- FALSE

    }else if (pophaz == "rescaled") {
      rescale <- TRUE
    }
    nalpha <- 0
  }

  if (rescale) {
    if (nTD >= 1) {
      name[6:(5 + (5 * nTD))] <- sapply(1:nTD,
                                        function(i, name.bs)
                                          rep({paste("qbs",
                                                     name.bs[i],
                                                     collapse = "")},
                                              5),
                                        name.bs = name.bs)
      value <- c(rep(1:5, nTD))
      name[6:(5 + (5 * nTD))] <- sapply(1:(5 * (nTD)),
                                        function(i, name, value) {
                                          paste(name[5 + i],
                                                "(", value[i], ")",
                                                collapse = "")
                                        },
                                        name = name,
                                        value = value)
      #name alpha
      if (nalpha == 1) {
        name[(5 + (5 * nTD) + 1):(5 + (5 * nTD) + nalpha)] <- "log(alpha)"
      }else{
        name[(5 + (5 * nTD) + 1):(5 + (5 * nTD) + nalpha)] <- paste0(
          'log(alpha.', levels(add.rmap), ")")
      }
      #name alpha
      if (nPH != 0) {
        name[(5 + (5 * nTD) + 1):(5 + (5 * nTD) + nPH)] <- name.bs[
          (nTD + 1):(nTD + nPH)]

        if (nalpha == 1) {
          name[(5 + 1 + (5 * nTD) + nPH):(nalpha + 5 + (5 * nTD) + nPH)] <-
            "log(alpha)"
        }else{
          name[
            (1 + 5 + (5 * nTD) + nPH):(nalpha + 5 + (5 * nTD) + nPH)
          ] <- paste0('log(alpha.', levels(add.rmap), ")")
        }
      }
    }
    else{
      if (nalpha == 1) {
        name[(5 + 1):(5 + nPH)] <- name.bs[(nTD + 1):(nTD + nPH)]
        name[(1 + 5 + nPH):(nalpha + 5 + nPH)] <- rep("alpha", nalpha)
      }else{
        name[(5 + 1):(5 + nPH)] <- name.bs[(nTD + 1):(nTD + nPH)]
        name[(5 + 1 + nPH):(nalpha + 5 + nPH)] <- paste0('alpha.',
                                                         levels(add.rmap))
      }
    }

    ##If no nalpha
  } else{
    if (nTD >= 1) {
      name[6:(5 + (5 * nTD))] <-
        sapply(1:nTD, function(i, name.bs)
          rep({
            paste("qbs", name.bs[i], collapse = "")
          }, 5), name.bs = name.bs)
      value <- c(rep(1:5, nTD))
      name[6:(5 + (5 * nTD))] <- sapply(1:(5 * (nTD)),
                                        function(i, name, value) {
                                          paste(name[5 + i],
                                                "(",
                                                value[i],
                                                ")",
                                                collapse = "")
                                        },
                                        name = name,
                                        value = value)

      if (nPH != 0) {
        name[(5 + (5 * nTD) + 1):(5 + (5 * nTD) + nPH)] <- name.bs[
          (nTD + 1):(nTD + nPH)]
      }
    }
    else{
      name[(5 + 1):(5 + nPH)] <- name.bs[(nTD + 1):(nTD + nPH)]
    }
  }
  covxx <- matrix(0, nrow(x), ncol(x) ^ 2)
  covxx <- matrix(
    sapply(1:ncol(x),
           function(i, xx)
             xx * xx[, i],
           xx = x),
    ncol = ncol(x) ^ 2)

  f <- function(theta0) {
    logLik <- 0
    rrBetaPH <- 1
    rrBetaTD <- 1
    rrTauTD <- exp(colSums(as.matrix(theta0[1:5] * t(P))))
    if (rescale) {
      levels(add.rmap) <- 1:nlevels(add.rmap)
      if (length(levels(add.rmap)) < 2) {
        alpha <- alpha0 <- (theta0[(5 + 5 * nTD + nPH + 1)])
      } else{
        no <- length(add.rmap)
        alpha0 <- (theta0[
          (5 + 5 * nTD + nPH + 1):(5 + 5 * nTD + nPH + nalpha)])
        Madd.rmap <- t(sapply(1:no, function(i, add.rmap) {
          ligne <- rep(0, nlevels(add.rmap))
          ligne[as.numeric(add.rmap[i])] <- 1
          return(ligne)
        }, add.rmap = add.rmap))
        alpha <- (Madd.rmap %*% alpha0)
      }
    } else{
      alpha <- 0
    }
    # rrBetaPH <- 1
    # rrBetaTD <- 1
    # rrTauTD <- exp(colSums(as.matrix(theta0[1:5] * t(P))))

    if (nPH != 0) {
      thetaPH <- theta0[(5 + 5 * nTD + 1):(5 + 5 * nTD + nPH)]
      rrBetaPH <- exp(colSums(as.matrix(
        thetaPH * t(x[, (1 + nTD):(nTD + nPH)]))))
    }
    if (nTD != 0) {
      rrBetaTD <- c(rep(0, nrowx))
      rrBetaTD <-
        exp(rowSums(
          sapply(1:5, function(i, theta, nTD, x, P)
            colSums(as.matrix(
              theta[(6 + (i - 1) * nTD):(5 + i * nTD)] *
                t(x[, 1:nTD] * P[, i]))), theta = theta0[1:(5 + 5 * nTD)],
            nTD = nTD, x = x, P = P)
        ))
    }
    chazardTDPH <- rrTauTD * rrBetaPH * rrBetaTD
    int.chazard <- matrix(0, nrowx, k)
    int.FDbase <- matrix(0, nrowx, 5 * k)
    int.SDbase <- matrix(0, 5 * nrowx, 5 * k)

    #Loop for integral calculation
    for (cst in 1:k) {
      int.chazard[, cst] <-
        IntGL(Int.chazard, t(matrix(
          c(boundmax[, cst], boundmin[, cst]),
          ncol = 2)),
          cst,
          cpti,
          cptj,
          theta0[1:(5 + 5 * nTD)],
          x,
          nTD,
          p, nghq = nghq)
    }

    if (rescale) {
      tlogLik <- colSums(as.matrix((-rrBetaPH * rowSums(int.chazard)) -
                                     exp(alpha) * ehazardInt +
                                     event[1:nrowx] *log(chazardTDPH +
                                                           exp(alpha) *
                                                           ehazard)))
    }
    else{
      tlogLik <- colSums(as.matrix(((-rrBetaPH * rowSums(int.chazard)) -
                                      exp(0) * ehazardInt +
                                      event[1:nrowx] *
                                      log(chazardTDPH + ehazard))))
    }
    logLik <- logLik + tlogLik
    return(logLik)
  }

  gradient <- function(theta0) {
    FD <-	FDbPH <- FDbase <- FDbTD <- 0
    FDalpha <- 0

    if (rescale) {
      levels(add.rmap) <- 1:nlevels(add.rmap)
      if (length(levels(add.rmap)) < 2) {
        alpha <- alpha0 <- (theta0[(5 + 5 * nTD + nPH + 1)])
      } else{
        no <- length(add.rmap)
        alpha0 <- (theta0[
          (5 + 5 * nTD + nPH + 1):(5 + 5 * nTD + nPH + nalpha)])
        Madd.rmap <- t(sapply(1:no, function(i, add.rmap) {
          ligne <- rep(0, nlevels(add.rmap))
          ligne[as.numeric(add.rmap[i])] <- 1
          return(ligne)
        }, add.rmap = add.rmap))
        alpha <- (Madd.rmap %*% alpha0)
      }
    } else{
      alpha <- 0
    }

    rrBetaPH <- 1
    rrBetaTD <- 1
    rrTauTD <- exp(colSums(as.matrix(theta0[1:5] * t(P))))

    if (nPH != 0) {
      thetaPH <- theta0[(5 + 5 * nTD + 1):(5 + (5 * nTD) + nPH)]
      rrBetaPH <- exp(colSums(as.matrix(
        thetaPH * t(x[, (1 + nTD):(nTD + nPH)]))))
    }

    if (nTD != 0) {
      rrBetaTD <- c(rep(0, nrowx))
      rrBetaTD <-
        exp(rowSums(
          sapply(1:5, function(i, theta, nTD, x, P)
            colSums(as.matrix(theta[(6 + (i - 1) * nTD):(5 + i * nTD)] *
                                t(x[, 1:nTD] * P[, i]))),
            theta = theta0[1:(5 + 5 * nTD)], nTD = nTD, x = x, P = P)
        ))
    }

    chazardTDPH <- rrTauTD * rrBetaPH * rrBetaTD
    int.chazard <- matrix(0, nrowx, k)
    int.FDbase <- matrix(0, nrowx, 5 * k)
    int.SDbase <- matrix(0, 5 * nrowx, 5 * k)

    #Loop for integral calculation
    for (cst in 1:k) {
      int.chazard[, cst] <-
        IntGL(Int.chazard, t(matrix(c(
          boundmax[, cst], boundmin[, cst]
        ), ncol = 2)), cst, cpti, cptj, theta0[1:(5 + 5 * nTD)], x,
        nTD, p, nghq = nghq)

      for (cpti in 1:5) {
        int.FDbase[, (cpti - 1) * k + cst] <- IntGL(
          Int.FDbase,
          t(matrix(c(boundmax[, cst],
                     boundmin[, cst]),
                   ncol = 2)),
          cst,
          cpti,
          cptj,
          theta0[1:(5 + 5 * nTD)],
          x,
          nTD,
          p, nghq = nghq)

        for (cptj in cpti:5) {
          int.SDbase[((cpti - 1) * nrowx + 1):(cpti * nrowx),
                     (cptj - 1) * k + cst] <- IntGL(Int.SDbase,
                                                    t(matrix(
                                                      c(boundmax[, cst],
                                                        boundmin[, cst]),
                                                      ncol = 2)),
                                                    cst,
                                                    cpti,
                                                    cptj,
                                                    theta0[1:(5 + 5 * nTD)],
                                                    x,
                                                    nTD,
                                                    p, nghq = nghq)

          int.SDbase[((cptj - 1) * nrowx + 1):(cptj * nrowx),
                     (cpti - 1) * k + cst] <- int.SDbase[(
                       (cpti - 1) * nrowx + 1):(cpti * nrowx),
                       (cptj - 1) * k + cst]
        }
      }
    }
    #First derivative for the baseline
    if (rescale) {
      FDbase <- sapply(1:5, function(i,
                                     event,
                                     int.FDbase,
                                     chazardTDPH,
                                     rrBetaPH,
                                     ehazard,
                                     nrowx,
                                     P)
        -sum(rrBetaPH * int.FDbase[, ((i - 1) * 3 + 1):(i * 3)]) +
          colSums(as.matrix(
            event[1:nrowx] * P[, i] *
              chazardTDPH / (chazardTDPH + exp(alpha) * ehazard))),
        event = event, int.FDbase = int.FDbase,
        chazardTDPH = chazardTDPH, rrBetaPH = rrBetaPH,
        ehazard = ehazard, nrowx = nrowx,
        P = P)
    } else{
      FDbase <-
        sapply(1:5, function(i,
                             event,
                             int.FDbase,
                             chazardTDPH,
                             rrBetaPH,
                             ehazard,
                             nrowx,
                             P)
          -sum(rrBetaPH * int.FDbase[, ((i - 1) * 3 + 1):(i * 3)]) +
            colSums(as.matrix(event[1:nrowx] * P[, i] *
                                chazardTDPH / (chazardTDPH + ehazard))),
          event = event, int.FDbase = int.FDbase,
          chazardTDPH = chazardTDPH, rrBetaPH = rrBetaPH,
          ehazard = ehazard, nrowx = nrowx, P = P)
    }

    if (nTD != 0) {
      # First derivative for the TD beta
      FDbTD <- lapply(1:nTD,
                      function(i, int.FDbase, x)
                        int.FDbase * x[, i]
                      ,int.FDbase = int.FDbase,
                      x = x)

      dum <- matrix(0, 5, nTD)
      for (j in 1:nTD) {
        dum[, j] <- sapply(1:5,
                           function(i, j, FDbTD, rrBetaPH)
                             -sum((rrBetaPH) * FDbTD[[j]][, (3 * (i - 1) + 1):(3 * i)]),
                           j = j,
                           FDbTD = FDbTD,
                           rrBetaPH = rrBetaPH)
      }

      if (rescale) {
        FDbTD <-
          c(t(dum) + matrix(
            sapply(1:5, function(i,
                                 event,
                                 chazardTDPH,
                                 x,
                                 nTD,
                                 ehazard,
                                 P)
              colSums(as.matrix(
                event * x[, 1:nTD] * P[, i] * chazardTDPH /
                  (chazardTDPH + exp(alpha) * ehazard))),
              ehazard = ehazard, event = event,
              chazardTDPH = chazardTDPH,
              x = x, nTD = nTD, P = P),
            nrow = nTD
          ))
      } else{
        FDbTD <- c(t(dum) + matrix(
          sapply(1:5, function(i,
                               event,
                               chazardTDPH,
                               x,
                               nTD,
                               ehazard,
                               P)
            colSums(as.matrix(event * x[, 1:nTD] * P[, i] *
                                chazardTDPH / (chazardTDPH + ehazard))),
            ehazard = ehazard,
            event = event,
            chazardTDPH = chazardTDPH,
            x = x,
            nTD = nTD, P = P),
          nrow = nTD
        ))
      }

      if (rescale) {
        if (length(levels(add.rmap)) < 2) {
          FDalpha <- colSums(as.matrix(
            exp(alpha) * (-ehazardInt) + event[1:nrowx] *
              c((exp(alpha)) * ehazard) / c(chazardTDPH +
                                              exp(alpha) * ehazard)))
        } else{
          Malpha <- t(alpha0 * t(Madd.rmap))
          FDalpha <- colSums(as.matrix(
            -exp(Malpha) * c(ehazardInt) + event[1:nrowx] *
              c(exp(Malpha) * ehazard) / c(chazardTDPH + exp(Malpha) *
                                             ehazard)))
        }
        FD <- c(FDbase , FDalpha, FDbTD)
      } else{
        FD <- c(FDbase, FDbTD)
      }

      #If there are some PH covariates associated
      if (nPH != 0) {
        if (rescale) {
          FDbPH <- colSums(as.matrix(
            (-x[, (1 + nTD):(nTD + nPH)] * rrBetaPH *
               rowSums(int.chazard) + event[1:nrowx] *
               x[, (1 + nTD):(nTD + nPH)] *
               (chazardTDPH / (chazardTDPH + (exp(alpha)) * ehazard)))))
          if (length(levels(add.rmap)) < 2) {
            FDalpha <- colSums(as.matrix(
              (exp(alpha)) * (-ehazardInt) + event[1:nrowx] *
                c((exp(alpha)) * ehazard) / c(chazardTDPH +
                                                (exp(alpha)) * ehazard)))
          } else{
            Malpha <- t(alpha0 * t(Madd.rmap))
            FDalpha <- colSums(as.matrix(
              -exp(Malpha) * c(ehazardInt) + event[1:nrowx] *
                c(exp(Malpha) * ehazard) / c(chazardTDPH +
                                               exp(Malpha) * ehazard)))
          }
          FD <- c(FDbase, FDbTD, FDbPH, FDalpha)

        } else{
          FDbPH <- colSums(as.matrix(
            -x[, (1 + nTD):(nTD + nPH)] *
              rrBetaPH * rowSums(int.chazard) +
              event * x[, (1 + nTD):(nTD + nPH)] *
              (chazardTDPH / (chazardTDPH + ehazard))))
          FD <- c(FDbase, FDbTD, FDbPH)
        }
      } else{
        if (rescale) {
          FD <- c(FDbase , FDalpha, FDbTD)
        } else{
          FD <- c(FDbase, FDbTD)
        }
      }
    }
    #If all the covariates are PH
    else{
      if (rescale) {
        if (length(levels(add.rmap)) < 2) {
          FDbPH <- colSums(as.matrix(
            -x[, (1 + nTD):(nTD + nPH)] * rrBetaPH *
              rowSums(int.chazard) + event * x[, (1 + nTD):(nTD + nPH)] *
              (chazardTDPH / (chazardTDPH + ehazard * exp(alpha)))))

          FDalpha <- colSums(as.matrix(
            -exp(alpha) * ehazardInt + event[1:nrowx] * exp(alpha) *
              c(ehazard) / c(chazardTDPH + exp(alpha) * ehazard)))

        } else{
          ##If there is more than j =1 alpha
          Malpha <- t(alpha0 * t(Madd.rmap))
          FDbPH <- colSums(as.matrix(-x[, (1 + nTD):(nTD + nPH)] * rrBetaPH *
                                       rowSums(int.chazard) +
                                       event * x * chazardTDPH /
                                       c(chazardTDPH +
                                           ehazard * exp(alpha))))

          FDalpha <- colSums(as.matrix(
            -exp(Malpha) * ehazardInt + event * exp(Malpha) *
              c(ehazard) / c(chazardTDPH + exp(alpha) * ehazard)))
        }
        FD <- c(FDbase, FDbPH, FDalpha)
      } else{
        FDbPH <- colSums(as.matrix(
          -x[, (1 + nTD):(nTD + nPH)] * rrBetaPH *
            rowSums(int.chazard) + event *
            x[, (1 + nTD):(nTD + nPH)] *
            (chazardTDPH / (chazardTDPH + ehazard))))
        FD <- c(FDbase, FDbPH)
      }
      FD
    }
  }

  #initialisation
  if (rescale) {

    if ((length(levels(add.rmap)) < 2) & (length(levels(add.rmap)) >= 1)) {
      nalphaLev <- 1
    } else{
      if (length(levels(add.rmap)) >= 2) {
        nalphaLev <- nlevels(add.rmap)
      }
    }

    theta0 <- theta0[1:(length(theta0) - nalphaLev)]

    nalpha <- 0
    rescale <- FALSE
    if (nTD != 0 & nPH != 0) {

      if (method == "L-BFGS-B") {
        theta0 <- optim(par = theta0,
                        fn = f,
                        gr = gradient,
                        hessian = TRUE,
                        lower = c(rep(-15, 5), rep(-15, (5) * nTD),
                                  rep((-15), nPH),
                                  rep(-15, nlevels(add.rmap))),
                        upper = c(rep(10, 5), rep(10, 5 * nTD), rep(10, nPH),
                                  rep(10, nlevels(add.rmap))),
                        method =  "L-BFGS-B",
                        control = list(REPORT = 1,
                                       maxit = 400,
                                       fnscale = -1,
                                       trace = trace))$par
      } else {
        theta0 <- optim(par = theta0,
                        fn = f,
                        hessian = TRUE,
                        method =  method,
                        control = list(REPORT = 1,
                                       maxit = 400,
                                       fnscale = -1,
                                       trace = trace))$par
      }

    } else{
      if (nTD == 0 & nPH != 0) {
        if (speedy) {
          max_cores <- parallel::detectCores()
          used_cores <- min(max_cores, (max_cores - 2))
          if (used_cores < 2) {
            stop("We didn't detect enough cores for speedy")
          }
          cl <- makeCluster(used_cores)
          setDefaultCluster(cl = cl)
          theta0 <- optimParallel(par = theta0,
                                  fn = f,
                                  gr = gradient,
                                  hessian = TRUE,
                                  method = method,
                                  control = list(REPORT = 1,
                                                 maxit = 400,
                                                 fnscale = -1,
                                                 trace = trace),
                                  parallel = list(loginfo = TRUE))
          setDefaultCluster(cl = NULL)
          stopCluster(cl)
        } else{
          if (method == "L-BFGS-B") {
            theta0 <- optim(par = theta0,
                            fn = f,
                            gr = gradient,
                            lower = c(rep(-15, 5), rep(-15, (5) * nTD),
                                      rep((-15), nPH), rep(-15, nalpha)),
                            upper = c(rep(10, 5), rep(10, 5 * nTD), rep(10, nPH),
                                      rep(15, nalpha)),
                            hessian = TRUE,
                            control = list(maxit = 400),
                            method = "L-BFGS-B")$par
          }else {
            theta0 <- optim(par = theta0,
                            fn = f,
                            hessian = TRUE,
                            control = list(maxit = 400),
                            method = method)$par
          }
        }
      } else {
        if (nTD != 0 & nPH == 0) {
          if (method == "L-BFGS-B") {
            theta0 <- optim(par = theta0,
                            fn = f,
                            gr = gradient,
                            hessian = TRUE,
                            lower = c(rep(-15, 5),
                                      rep(-15 * nTD),
                                      rep(-15, nalpha)),
                            upper = c(rep(15, 5),
                                      rep(15, 5 * nTD),
                                      rep(15, nalpha)),
                            method = "L-BFGS-B",
                            control = list(REPORT = 1,
                                           maxit = 400,
                                           fnscale = -1,
                                           trace = trace))$par
          }else {
            theta0 <- optim(par = theta0,
                            fn = f,
                            hessian = TRUE,
                            method = method,
                            control = list(REPORT = 1,
                                           maxit = 400,
                                           fnscale = -1,
                                           trace = trace))$par
          }

        }
      }
    }
  } else {
    if (nTD != 0 & nPH != 0) {
      if (method == "L-BFGS-B") {
        theta0 <- optim(par = theta0,
                        fn = f,
                        gr = gradient,
                        hessian = TRUE,
                        lower = c(rep(-15, 5),
                                  rep(-15, (5) * nTD),
                                  rep((-15), nPH)),
                        upper = c(rep(10, 5),
                                  rep(10, 5 * nTD),
                                  rep(10, nPH)),
                        method =  "L-BFGS-B",
                        control = list(REPORT = 1,
                                       maxit = 400,
                                       fnscale = -1,
                                       trace = trace))$par
      }else {
        theta0 <- optim(par = theta0,
                        fn = f,
                        hessian = TRUE,
                        method = method,
                        control = list(REPORT = 1,
                                       maxit = 400,
                                       fnscale = -1,
                                       trace = trace))$par
      }

    }
    else{

      if (nTD == 0 & nPH != 0) {
        if (method == "L-BFGS-B") {
          theta0 <- optim(par = theta0,
                          fn = f,
                          gr = gradient,
                          hessian = TRUE,
                          method = "L-BFGS-B",
                          lower = c(rep(-15, 5), rep((-15), nPH)),
                          upper = c(rep(15, 5), rep(15, nPH)),
                          control = list(REPORT = 1,
                                         maxit = 500,
                                         fnscale = -1,
                                         trace = trace))$par
        }else {
          theta0 <- optim(par = theta0,
                          fn = f,
                          hessian = TRUE,
                          method = method,
                          control = list(REPORT = 1,
                                         maxit = 500,
                                         fnscale = -1,
                                         trace = trace))$par
        }

      }else{
        if (nTD != 0 & nPH == 0) {
          if (method == "L-BFGS-B") {
            theta0 <- optim(par = theta0,
                            fn = f,
                            gr = gradient,
                            hessian = TRUE,
                            lower = c(rep(-15), rep(-15 * nTD)),
                            upper = c(rep(15, 10), rep(15, 5 * nTD)),
                            method =  "L-BFGS-B",
                            control = list(REPORT = 1,
                                           maxit = 500,
                                           fnscale = -1,
                                           trace = trace))$par
          }else {
            theta0 <- optim(par = theta0,
                            fn = f,
                            hessian = TRUE,
                            method = method,
                            control = list(REPORT = 1,
                                           maxit = 500,
                                           fnscale = -1,
                                           trace = trace))$par
          }

        }
      }
    }
  }



  if (rescale == FALSE) {
    if (nTD != 0 & nPH != 0) {
      if (method == "L-BFGS-B") {
        res <- optim(par = theta0,
                     fn = f,
                     gr = gradient,
                     hessian = TRUE,
                     method =  "L-BFGS-B",
                     control = list(REPORT = 1,
                                    maxit = 5000,
                                    fnscale = -1,
                                    trace = trace))
      }else {
        res <- optim(par = theta0,
                     fn = f,
                     hessian = TRUE,
                     method =  method,
                     control = list(REPORT = 1,
                                    maxit = 5000,
                                    fnscale = -1,
                                    trace = trace))
      }

    }else{
      if (nTD == 0 & nPH != 0) {
        if (speedy) {
          max_cores <- parallel::detectCores()
          used_cores <- min(max_cores, (max_cores - 2))
          if (used_cores < 2) {
            stop("We didn't detect enough cores for speedy")
          }
          cl <- makeCluster(used_cores)
          setDefaultCluster(cl = cl)
          res <- optimParallel(par = theta0,
                               fn = f,
                               gr = gradient,
                               hessian = TRUE,
                               method = method,
                               control = list(REPORT = 1,
                                              maxit = 5000,
                                              fnscale = -1,
                                              trace = trace),
                               parallel = list(loginfo = TRUE))
          setDefaultCluster(cl = NULL)
          stopCluster(cl)
        }
        else{
          if (method == "L-BFGS-B") {
            res <- optim(par = theta0,
                         fn = f,
                         gr = gradient,
                         hessian = TRUE,
                         method = "L-BFGS-B",
                         control = list(REPORT = 1,
                                        maxit = 5000,
                                        fnscale = -1,
                                        trace = trace))
          } else {
            res <- optim(par = theta0,
                         fn = f,
                         hessian = TRUE,
                         method = method,
                         control = list(REPORT = 1,
                                        maxit = 5000,
                                        fnscale = -1,
                                        trace = trace))
          }
        }
      }else {
        if (nTD != 0 & nPH == 0) {
          if (method == "L-BFGS-B") {
            res <- optim( par = theta0,
                          fn = f,
                          gr = gradient,
                          hessian = TRUE,
                          method =  "L-BFGS-B",
                          control = list(REPORT = 1,
                                         maxit = 5000,
                                         fnscale = -1,
                                         trace = trace))
          } else {
            res <- optim( par = theta0,
                          fn = f,
                          hessian = TRUE,
                          method =  method,
                          control = list(REPORT = 1,
                                         maxit = 5000,
                                         fnscale = -1,
                                         trace = trace))
          }
        }
      }
    }
  }

  if ((length(levels(add.rmap)) < 2) & length(levels(add.rmap)) >= 1) {
    nalpha <- 1
    rescale <- TRUE
  } else{
    if (length(levels(add.rmap)) >= 2) {
      rescale <- TRUE
      nalpha <- nlevels(add.rmap)
    } else{
      if (length(levels(add.rmap)) < 1) {
        rescale <- FALSE
        nalpha <- 0
      }
    }
  }

  if (rescale) {
    theta0 <- c(res$par, rep(0.1, nalpha))
    if (nTD != 0 & nPH != 0) {
      if (method == "L-BFGS-B") {
        res <- optim(par = theta0,
                     fn = f,
                     gr = gradient,
                     hessian = TRUE,
                     method =  "L-BFGS-B",
                     control = list(REPORT = 1,
                                    maxit = 5000,
                                    fnscale = -1,
                                    trace = trace))
      }else {
        res <- optim(par = theta0,
                     fn = f,
                     hessian = TRUE,
                     method =  method,
                     control = list(REPORT = 1,
                                    maxit = 5000,
                                    fnscale = -1,
                                    trace = trace))
      }
    }
    else{
      if (nTD == 0 & nPH != 0) {
        if (speedy) {
          max_cores <- parallel::detectCores()
          used_cores <- min(max_cores, (max_cores - 2))
          if (used_cores < 2) {
            stop("We didn't detect enough cores for speedy")
          }
          cl <- parallel::makeCluster(used_cores)
          setDefaultCluster(cl = cl)
          res <- optimParallel(par = theta0,
                               fn = f,
                               gr = gradient,
                               hessian = TRUE,
                               method = "method",
                               control = list(REPORT = 1,
                                              maxit = 5000,
                                              fnscale = -1,
                                              trace = trace),
                               parallel = list(loginfo = TRUE))
          setDefaultCluster(cl = NULL)
          stopCluster(cl)
        }
        else{
          if (method == "L-BFGS-B") {
            res <- optim(par = theta0,
                         fn = f,
                         gr = gradient,
                         hessian = TRUE,
                         method = "L-BFGS-B",
                         control = list(REPORT = 1,
                                        maxit = 5000,
                                        fnscale = -1,
                                        trace = trace))
          }else {
            res <- optim(par = theta0,
                         fn = f,
                         hessian = TRUE,
                         method = method,
                         control = list(REPORT = 1,
                                        maxit = 5000,
                                        fnscale = -1,
                                        trace = trace))
          }
        }
      } else{
        if (nTD != 0 & nPH == 0) {
          if (method == "L-BFGS-B") {
            res <- optim(par = theta0,
                         fn = f,
                         gr = gradient,
                         hessian = TRUE,
                         method =  "L-BFGS-B",
                         control = list(REPORT = 1,
                                        maxit = 5000,
                                        fnscale = -1,
                                        trace = FALSE))
          }else {
            res <- optim(par = theta0,
                         fn = f,
                         hessian = TRUE,
                         method =  method,
                         control = list(REPORT = 1,
                                        maxit = 5000,
                                        fnscale = -1,
                                        trace = FALSE))
          }
        }
      }
    }
  }





  logLik <- res$value
  theta0 <- res$par
  FD <- -gradient(theta0)
  SD <- -numDeriv::hessian(f, theta0)
  iter <- res$counts[2]
  convergence <- res$convergence
  message <- res$message
  names(theta0) <- name

  return(list(coefficients = theta0,
              varcov = try(solve(SD), TRUE),
              std_err = try(sqrt(diag(solve(SD))), TRUE),
              loglik = logLik,
              iterations = iter,
              intervalles = int,
              convergence = convergence,
              message = message,
              nTD = nTD,
              nPH = nPH,
              nalpha = nalpha))
}








