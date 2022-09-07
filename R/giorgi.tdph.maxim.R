

giorgi.tdph.maxim <- function(x, theta0, nTD, nPH, event, ehazard, ehazardInt,
                              P, p, k, nrowx, IntGL, nghq = nghq, cpti, cptj, cst, boundmin,
                              boundmax, Int.chazard, Int.FDbase, Int.SDbase,
                              int, control, speedy =  FALSE) {
  control.iter.max <- control$iter.max
  control.eps <- control$eps
  name.bs <- dimnames(x)[[2]]
  name <- c(rep(0, 5 + 5 * (nTD) + nPH))
  name[1:5] <- c("qbs base ( 1 )", "qbs base ( 2 )", "qbs base ( 3 )",
                 "qbs base ( 4 )", "qbs base ( 5 )")

  if (nTD >= 1) {
    name[6:(5 + (5 * nTD))] <- sapply(1:nTD,
                                      function(i, name.bs)
                                        rep(
                                          {paste("qbs",
                                                 name.bs[i],
                                                 collapse = "") },
                                          5),
                                      name.bs = name.bs)
    value <- c(rep(1:5, nTD))
    name[6:(5 + (5 * nTD))] <- sapply(1:(5 * (nTD)),
                                      function(i, name, value) {
                                        paste(name[5 + i],
                                              "(", value[i], ")",
                                              collapse = "")},
                                      name = name, value = value)
    if (nPH != 0)
      name[(5 + (5 * nTD) + 1):(5 + (5 * nTD) + nPH)] <-
      name.bs[(nTD + 1):(nTD + nPH)]
  }
  else
    name[(5 + 1):(5 + nPH)] <- name.bs[(nTD + 1):(nTD + nPH)]

  covxx <- matrix(0, nrow(x), ncol(x) ^ 2)
  covxx <- matrix(sapply(1:ncol(x),
                         function(i, xx)
                           xx * xx[, i], xx = x),
                  ncol = ncol(x) ^ 2)
  iter <- 0
  FD <- 1
  SD <- NULL
  #Maximlisation loop
  while (sum(abs(FD)) > (control.eps))
  {
    if (iter > control.iter.max)
      stop(paste("Ran out of iterations", control.iter.max,
                 "and did not converge." ))
    iter <- iter + 1
    #Calculation of the disease-related mortality hazard function with the
    #TD parameters, and with the PH parameters if the are some
    rrBetaPH <- 1
    rrBetaTD <- 1
    rrTauTD <- exp(colSums(theta0[1:5] * t(P)))
    if (nPH != 0) {
      thetaPH <- theta0[(5 + 5 * nTD + 1):(5 + 5 * nTD + nPH)]
      rrBetaPH <-
        exp(colSums(as.matrix((thetaPH * t(x[, (1 + nTD):(nTD + nPH)])))))
    }
    if (nTD != 0) {
      rrBetaTD <- c(rep(0, nrowx))
      rrBetaTD <-
        exp(rowSums(
          sapply(1:5, function(i, theta, nTD, x, P)
            colSums(as.matrix((theta[(6 + (i - 1) * nTD):(5 + i * nTD)] *
                                 t(x[, 1:nTD] * P[, i])))),
            theta = theta0[1:(5 + 5 * nTD)],
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
        IntGL(Int.chazard, t(matrix(c(
          boundmax[, cst], boundmin[, cst]
        ), ncol = 2)), cst, cpti, cptj, theta0[1:(5 + 5 * nTD)], x, nTD, p, nghq = nghq)
      for (cpti in 1:5) {
        int.FDbase[, (cpti - 1) * k + cst] <-
          IntGL(Int.FDbase,
                t(matrix(c(boundmax[, cst], boundmin[, cst]), ncol = 2)),
                cst, cpti, cptj,
                theta0[1:(5 + 5 * nTD)],
                x,
                nTD,
                p, nghq)
        for (cptj in cpti:5) {
          int.SDbase[((cpti - 1) * nrowx + 1):(cpti * nrowx),
                     (cptj - 1) * k + cst] <-
            IntGL(Int.SDbase,
                  t(matrix(c(boundmax[, cst], boundmin[, cst] ), ncol = 2)),
                  cst, cpti, cptj,
                  theta0[1:(5 + 5 * nTD)],
                  x,
                  nTD,
                  p, nghq = nghq)
          int.SDbase[((cptj - 1) * nrowx + 1):(cptj * nrowx),
                     (cpti - 1) * k + cst] <-
            int.SDbase[((cpti - 1) * nrowx + 1):(cpti * nrowx),
                       (cptj - 1) * k + cst]
        }
      }
    }
    #First derivative for the baseline
    FDbase <-
      sapply(1:5, function(i,
                           event,
                           int.FDbase,
                           chazardTDPH,
                           rrBetaPH,
                           ehazard,
                           nrowx,
                           P)
        - sum(rrBetaPH * int.FDbase[, ((i - 1) * 3 + 1):(i * 3)]) +
          colSums(as.matrix(
            (event[1:nrowx] * P[, i] * chazardTDPH /
               (chazardTDPH + ehazard)))),
        event = event, int.FDbase = int.FDbase, chazardTDPH = chazardTDPH,
        rrBetaPH = rrBetaPH, ehazard = ehazard,
        nrowx = nrowx, P = P)

    #Second derivative for the baseline
    SDbase <- sapply(1:5, function(i, rrBetaPH, int.SDbase)
        (rrBetaPH) * rowSums(int.SDbase[, ((i - 1) * 3 + 1):(i * 3)]),
        rrBetaPH = rrBetaPH, int.SDbase = int.SDbase)
    SDbase <- sapply(1:5,
                     function(i, SDbase, nrowx)
                       colSums(
                         as.matrix(
                           (SDbase[
                             ((i - 1) * nrowx + 1):((i - 1) * nrowx + nrowx),
                             ]))),
                     SDbase = SDbase, nrowx = nrowx)
    covpp <- matrix(0, nrow(x), ncol(P) ^ 2)
    covpp <- matrix(sapply(1:ncol(P), function(i, PP)
      PP * PP[, i], PP = P), ncol = ncol(P) ^ 2)
    SDbase <- SDbase - colSums(event[1:nrowx] * covpp * chazardTDPH *
                               ehazard / (chazardTDPH + ehazard) ^ 2)

    #ATTENTION: the matrix of the second derivative is ordered in this way:
    # - components of the baseline
    # - components of the TD coavriates: columns for the first elementary
    # spline for the nTD covariates
    #   following by the columns for the second elementary spline for the nTD
    #   covariates, and so on
    #   for the 5 elementary splines
    # - components of the PH covariates if necessary

    #If there are some TD covariates
    if (nTD != 0) {
      #First derivative for the TD beta
      FDbTD <- lapply(1:nTD, function(i, int.FDbase, x)
          int.FDbase * x[, i], int.FDbase = int.FDbase, x = x)
      dum <- matrix(0, 5, nTD)
      for (j in 1:nTD) {
        dum[, j] <- sapply(1:5,
                           function(i, j, FDbTD, rrBetaPH)
                             - sum(
                               (rrBetaPH) *
                                 FDbTD[[j]][,(3 * (i - 1) + 1):(3 * i)]),
                           j = j,
                           FDbTD = FDbTD,
                           rrBetaPH = rrBetaPH)
      }

      FDbTD <- c(t(dum) + matrix(
        sapply(1:5,
               function(i,
                        event,
                        chazardTDPH,
                        x,
                        nTD,
                        ehazard,
                        P)
            colSums(as.matrix((event * x[, 1:nTD] * P[, i] * chazardTDPH /
                                 (chazardTDPH + ehazard)))),
            ehazard = ehazard,
            event = event,
            chazardTDPH = chazardTDPH,
            x = x,
            nTD = nTD,
            P = P),
          nrow = nTD
        ))
      FD <- c(FDbase, FDbTD)

      #Second partial derivative
      SDP <-
        matrix(
          sapply(1:nTD, function(i, int.SDbase, x, nTD)
            int.SDbase * x[, i], int.SDbase = int.SDbase, x = x, nTD = nTD),
          ncol = ncol(int.SDbase) * nTD
        )
      dummy1 <-
        sapply(1:(5 * nTD), function(i, SDP, rrBetaPH)
          (rrBetaPH) * rowSums(SDP[, ((i - 1) * 3 + 1):(i * 3)]), SDP = SDP,
          rrBetaPH = rrBetaPH)
      dummy1 <-
        lapply(1:5, function(i, dummy1, nrowx, nTD)
          matrix(colSums(as.matrix(
            (dummy1[((i - 1) * nrowx + 1):((i - 1) * nrowx + nrowx), ]))),
            ncol = nTD), dummy1 = dummy1, nrowx = nrowx, nTD = nTD)
      dummy2 <- sapply(1:nTD, function(i,
                               event,
                               chazardTDPH,
                               covpp,
                               ehazard,
                               nrowx,
                               x)
          (-colSums(as.matrix((event[1:nrowx] * x[, i] * covpp * chazardTDPH *
                                 ehazard / (chazardTDPH + ehazard) ^ 2)))),
          event = event, chazardTDPH = chazardTDPH, x = x, nrowx = nrowx,
          ehazard = ehazard, covpp = covpp)

      SDP <- dummy2 + do.call("rbind", dummy1)
      SDPTD <-
        matrix(sapply(1:5, function(i, SDP)
          (SDP[((i - 1) * 5 + 1):(i * 5), ]), SDP = SDP), 5, 5 * nTD)

      #Second derivative for the TD beta
      covxxTD <- matrix(0, nrow(x), (nTD ^ 2))
      covxxTD <-
        matrix(sapply(1:(nTD), function(i, xx, nTD)
          xx[, 1:nTD] * xx[, i], xx = x, nTD = nTD),
          ncol = (nTD ^ 2))
      SDbTD <-
        matrix(
          sapply(1:(nTD ^ 2), function(i, int.SDbase, covxxTD)
            int.SDbase * covxxTD[, i], int.SDbase = int.SDbase,
            covxxTD = covxxTD),
          ncol = ncol(int.SDbase) * (nTD ^ 2)
        )
      n2 <- sapply(1:(5 * nTD ^ 2), function(i, SDbTD, rrBetaPH)
          (rrBetaPH) * rowSums(SDbTD[, ((i - 1) * 3 + 1):(i * 3)]),
          SDbTD = SDbTD,
          rrBetaPH = rrBetaPH)

      n1 <- lapply(1:5, function(i, n2, nrowx, nTD)
        matrix(colSums(as.matrix(
          (n2[((i - 1) * nrowx + 1):((i - 1) * nrowx + nrowx), ]))),
          ncol = (nTD ^ 2)), n2 = n2, nrowx = nrowx, nTD = nTD)

      n3 <- sapply(1:ncol(covpp), function(i,
                                       event,
                                       chazardTDPH,
                                       nrowx,
                                       covxxTD,
                                       covpp,
                                       ehazard)
          (-colSums(as.matrix((event[1:nrowx] * covxxTD * covpp[, i] *
                                 chazardTDPH * ehazard /
                                 (chazardTDPH + ehazard) ^ 2)))),
          event = event, chazardTDPH = chazardTDPH, ehazard = ehazard,
          nrowx = nrowx, covxxTD = covxxTD, covpp = covpp)

      SDbTD <- t(n3) + c(do.call("rbind", n1))
      if (nTD == 1) {
        SDbTD <- matrix(SDbTD, 5, 5)
      }
      else{
        SDbTD <-
          do.call("rbind",
                  lapply(1:25, function(i, SDbTD, nTD)
                    matrix(c(SDbTD[i, ]), ncol = nTD, nrow = nTD),
                    SDbTD = SDbTD, nTD = nTD))
        SDbTD <-
          matrix(sapply(1:5, function(i, SDbTD, nTD)
            (SDbTD[((i - 1) * (5 * nTD) + 1):(i * (5 * nTD)), ]),
            SDbTD = SDbTD, nTD = nTD),
            5 * nTD,
            5 * nTD)
        NULL
      }
      SD <- rbind(cbind(SDbase, SDPTD), cbind(t(SDPTD), SDbTD))

      #If there are some PH covariates associated
      if (nPH != 0) {
        FDbPH <- colSums(as.matrix((-x[, (1 + nTD):(nTD + nPH)] * rrBetaPH *
                                      rowSums(int.chazard) +
                                      event * x[, (1 + nTD):(nTD + nPH)] *
                                      (chazardTDPH /
                                         (chazardTDPH + ehazard)))))

        FD <- c(FD, FDbPH)

        SDPPH <- sapply(1:5, function(i, int.FDbase)
            rowSums(int.FDbase[, (3 * (i - 1) + 1):(3 * i)]),
            int.FDbase = int.FDbase)
        SDPPH <- sapply(1:nPH, function(i, SDPPH, rrBetaPH, x, nTD)
            - colSums(as.matrix(((-rrBetaPH) * x[, nTD + i] * SDPPH))),
            SDPPH = SDPPH, rrBetaPH = rrBetaPH, x = x, nTD = nTD)
        SDPPH1 <-
          sapply(1:nPH, function(i,
                                 event,
                                 chazardTDPH,
                                 nrowx,
                                 x,
                                 nTD,
                                 P,
                                 ehazard)
            - colSums(as.matrix((event[1:nrowx] * x[, nTD + i] * P *
                                   chazardTDPH *ehazard /
                                   (chazardTDPH +  ehazard) ^ 2))),
            event = event, chazardTDPH = chazardTDPH, nrowx = nrowx,
            x = x, nTD = nTD, P = P, ehazard = ehazard)
        SDPPH <- SDPPH + SDPPH1

        covxxTDPH <- matrix(0, nrowx, (nTD * nPH))
        covxxTDPH <- matrix(sapply(1:(nTD), function(i, xx, nTD, nPH)
            xx[, (nTD + 1):(nTD + nPH)] * xx[, i], xx = x,
            nTD = nTD, nPH = nPH),
            ncol = (nTD * nPH))
        SDPbTDPH1 <- matrix(t(sapply(1:ncol(covxxTDPH),
                                     function(i,
                                              event,
                                              chazardTDPH,
                                              nrowx,
                                              covxxTDPH,
                                              ehazard,
                                              P)
              - colSums(as.matrix((event[1:nrowx] * covxxTDPH[, i] * P *
                                     chazardTDPH * ehazard /
                                     (chazardTDPH + ehazard) ^ 2))),
              event = event, chazardTDPH = chazardTDPH, nrowx = nrowx,
              covxxTDPH = covxxTDPH, ehazard = ehazard, P = P)
          ), nPH, 5 * nTD)
        int.FDbTD <- matrix(
            sapply(1:nTD, function(i, int.FDbase, x)
              int.FDbase * x[, i], int.FDbase = int.FDbase, x = x),
            ncol = nTD * ncol(int.FDbase)
          )
        int.FDbTD <- matrix(
          sapply(1:(5 * nTD),
                 function(i, int.FDbTD, rrBetaPH, x, nTD, nPH)
              (-colSums(as.matrix(
                (-x[, (nTD + 1):(nTD + nPH)] * (rrBetaPH) *
                   rowSums(int.FDbTD[, ((i - 1) * 3 + 1):(i * 3)]))))),
              int.FDbTD = int.FDbTD, rrBetaPH = rrBetaPH, x = x,
              nTD = nTD, nPH = nPH),
          ncol = nTD)

        SDPbTDPH <- SDPbTDPH1 +
          c(do.call("cbind",
                    lapply(1:5, function(i, int.FDbTD, nPH)
                      int.FDbTD[((i - 1) * nPH + 1):((i - 1) * nPH + nPH), ],
                      int.FDbTD = int.FDbTD, nPH)))

        covxxPH <- matrix(0, nrowx, (nPH ^ 2))
        covxxPH <- matrix(sapply(1:(nPH), function(i, xx, nTD, nPH)
          xx[, (nTD + 1):(nTD + nPH)] * xx[, nTD + i],
          xx = x,
          nTD = nTD,
          nPH = nPH),
          ncol = (nPH ^ 2))
        SDbPH <- (-colSums(as.matrix(
          (-covxxPH * rrBetaPH * rowSums(int.chazard) +
             (event * covxxPH * (chazardTDPH * ehazard) /
                (chazardTDPH + ehazard) ^ 2)))))

        SDbPH <- matrix(SDbPH, nPH, nPH)

        SD <- rbind(cbind(SDbase, SDPTD, SDPPH),
                    (cbind(t(SDPTD), SDbTD, t(SDPbTDPH))),
                    (cbind(t(SDPPH), SDPbTDPH, SDbPH)))
        NULL
      }
    }
    #If all the covariates are PH
    else{
      FDbPH <- colSums(-x[, (1 + nTD):(nTD + nPH)] * rrBetaPH *
                  rowSums(int.chazard) +
                  event * x[, (1 + nTD):(nTD + nPH)] *
                  (chazardTDPH / (chazardTDPH + ehazard)))
      FD <- c(FDbase, FDbPH)

      SDPPH <-
        sapply(1:5, function(i, int.FDbase)
          rowSums(int.FDbase[, (3 * (i - 1) + 1):(3 * i)]),
          int.FDbase = int.FDbase)
      SDPPH <-
        sapply(1:nPH, function(i, SDPPH, rrBetaPH, x, nTD)
          - colSums(as.matrix(((-rrBetaPH) * x[, nTD + i] * SDPPH))),
          SDPPH = SDPPH, rrBetaPH = rrBetaPH, x = x, nTD = nTD)
      SDPPH1 <-
        sapply(1:nPH, function(i,
                               event,
                               chazardTDPH,
                               nrowx,
                               x,
                               nTD,
                               P,
                               ehazard)
          - colSums(as.matrix((event[1:nrowx] * x[, nTD + i] * P *
                                 chazardTDPH * ehazard /
                                 (chazardTDPH + ehazard) ^ 2))),
          event = event, chazardTDPH = chazardTDPH, nrowx = nrowx, x = x,
          nTD = nTD, P = P, ehazard = ehazard)
      SDPPH <- SDPPH + SDPPH1

      covxxPH <- matrix(0, nrowx, (nPH ^ 2))
      covxxPH <-
        matrix(sapply(1:(nPH), function(i, xx, nTD, nPH)
          xx[, (nTD + 1):(nTD + nPH)] * xx[, nTD + i], xx = x,
          nTD = nTD, nPH = nPH),
          ncol = (nPH ^ 2))
      SDbPH <-
        (-colSums(as.matrix(
          (-covxxPH * rrBetaPH * rowSums(int.chazard) +
             (event * covxxPH * (chazardTDPH * ehazard) /
                (chazardTDPH + ehazard) ^ 2)))))
      SDbPH <- matrix(SDbPH, nPH, nPH)

      SD <- rbind(cbind(SDbase, SDPPH), (cbind(t(SDPPH), SDbPH)))
      NULL
    }
    diff <- try(solve(qr(SD), FD))
    if (!is.numeric(diff))
  stop("Matrix not definite positive. Check for colinearity in the data set.")

    theta0 <- diff + theta0
    NULL
  }
  logLik <- colSums(as.matrix(((-rrBetaPH * rowSums(int.chazard)) -
                                        exp(0) * ehazardInt +
                                        (event[1:nrowx] * log(chazardTDPH +
                                                                ehazard)))))

  if (nTD != 0) {
    v <- c(rep(0, (5 * nTD)))
    for (i in 1:nTD) {
      v[((i - 1) * 5 + 1):((i - 1) * 5 + 5)] <- sapply(1:5,
                                                       function(j,
                                                                i,
                                                                theta0,
                                                                nTD)
                                                         theta0[
                                                           (6 + (j - 1) * nTD) +
                                                             (i - 1)],
                                                       i = i,
                                                       theta0 = theta0,
                                                       nTD)
    }
    theta0[6:(5 + (5 * nTD))] <- v
  }
  names(theta0) <- name
  list(
    coefficients = theta0,
    varcov = try(solve(SD), TRUE),
    std_err = try(sqrt(diag(solve(SD))), TRUE),
    loglik = logLik,
    iterations = iter,
    intervalles = int,
    nPH = nPH,
    nTD = nTD
  )
}
