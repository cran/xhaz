
exphaz_years <- function(ageDiag,
                         time,
                         data,
                         rmap,
                         ratetable,
                         ratedata = NULL,
                         add.rmap = NULL,
                         varlist,
                         temp01,
                         scale,
                         pophaz,
                         only_ehazard) {

  if (missing(ratetable)) {
    ageDC <- ageDiag + floor(time)
    ageDC <- (ageDC >= nrow(ratedata)) * (nrow(ratedata) - 1) +
      (ageDC < nrow(ratedata)) * ageDC
    ehazard <- (data[, rmap$sex] == "male") *
      ratedata[(ageDC + 1), 1] + (data[, rmap$sex] == "female") *
      ratedata[(ageDC + 1), 2]
    ratedataInt <- apply(ratedata, 2, cumsum)

    if (!is.null(add.rmap)) {
      ehazardInt_u <- (data[, rmap$sex] == "male") * (ageDC >= 1) *
        ratedataInt[ageDC, 1] +
        ratedata[(ageDC + 1), 1] * (ageDC - trunc(ageDC)) +
        (data[, rmap$sex] == "female") * (ageDC >= 1) *
        (ratedataInt[ageDC, 2]) +
        ratedata[ageDC + 1, 2] * (ageDC - trunc(ageDC))

      ehazardInt_l <-
        (data[, rmap$sex] == "male") * (ageDiag >= 1) *
        (ratedataInt[ageDiag, 1]) +
        ratedata[ageDiag + 1, 1] * (ageDiag - trunc(ageDiag))  +
        (data[, rmap$sex] == "female") * (ageDiag >= 1) *
        (ratedataInt[ageDiag, 2]) +
        ratedata[ageDiag + 1, 2] * (ageDiag - trunc(ageDiag))

      ehazardInt <- ehazardInt_u - ehazardInt_l
    }

    return(list(ehazard = ehazard,
                ehazardInt = ehazardInt))


  } else{
    if (!missing(ratedata)) {
      stop("Don't provide ratedata if ratetable is available!")
    }

    ageDC <- ageDiag + time
    if (max(ageDC) > 115) {
      stop(
        "Please check the scale used for time of follow-up: it must be in year; also check the age of diagnosis for older patients.\n"
      )
    }


    if (!inherits(data[, rmap$year], "Date")) {
      stop("Year must be in a Date class")
    }
    year <- as.integer(format(((data[, rmap$year]) +
                                 as.difftime(time * scale,
                                             units = "days")
    ), "%Y"))

    RT <- attr(ratetable, which = "dimnames")
    RTyear <- which(varlist == 'year')
    RTageDC <- which(varlist == 'age')

    if (max(year) > as.numeric(max(RT[[RTyear]]))) {
      idyearmax <- which(c(year > as.numeric(max(RT[[RTyear]]))))
      year[c(year > as.numeric(max(RT[[RTyear]])))] <- as.numeric(max(RT[[which(varlist == 'year')]]))
      cat("\nmax year+time in the data:\n", as.numeric(max(RT[[RTyear]])))
      cat("\nmax year+time considered:\n", as.numeric(max(RT[[RTyear]])))
    }


    if (max(ageDC) > as.numeric(max(RT[[RTageDC]]))) {
      idageDCmax <- which(c(ageDC > as.numeric(max(RT[[RTageDC]]))))

      ageDC[c(ageDC > as.numeric(max(RT[[RTageDC]])))] <- as.numeric(max(RT[[which(varlist == 'age')]]))
    cat("\nmax age+time in the data:\n", as.numeric(max(RT[[RTageDC]])))
    cat("\nmax age+time considered:\n", as.numeric(max(RT[[RTageDC]])))
      }
    nb.age.dc <- length(RT[[which(varlist == 'age')]])

    newdata01 <- data[, as.vector(unlist(rmap))[temp01]]
    names(newdata01) <- varlist
    newdata01$age <- trunc(ageDC)
    newdata01$year <- year

    Indic01  <- sapply(1:ncol(newdata01), function(i)
      apply(outer(newdata01[, temp01[i]], RT[[temp01[i]]], "=="), 1,
            function(x)
              which(x)))
    Indic01 <- Indic01[, temp01]
    colnames(Indic01) <- names(newdata01)

if (is.list((Indic01)[1,])) {
  ehazard <- mapply(function(i) {
    return( ratetable[matrix(unlist(as.data.frame(Indic01)[i,]),nrow = 1)])
  }, 1:length(ageDC)) * scale

}else if (is.integer((Indic01)[1,])) {
  ehazard <- mapply(function(i) {
    return(ratetable[matrix(Indic01[i, ], nrow = 1)])
  }, 1:length(ageDC)) * scale
}


    if (only_ehazard) {
      dateDiag <- data[, rmap$year]
      return(list(ehazard = ehazard,
                  dateDiag = dateDiag))
    } else{
      newdata02 <- data[, as.vector(unlist(rmap))[temp01]]#table matching life table and data
      names(newdata02) <- varlist
      newdata02$age <- trunc(ageDiag)
      newdata02$year <- as.integer(format(data[, rmap$year], '%Y'))



      Indic02 <-sapply(1:length(temp01), function(i) {
        apply(outer(newdata02[, temp01[i]],
                    RT[[temp01[i]]], "=="), 1,
              function(x) {
                which(x)
                })})

      Indic02 <- Indic02[, temp01]
      colnames(Indic02) <- names(newdata02)
      dateDC <- data[, rmap$year] + as.difftime(time * scale, units = "days")
      diffAge <- ceiling(data[, rmap$age]) - data[, rmap$age]
      dateAnniv <- data[, rmap$year] +
        as.difftime(diffAge * scale, units = "days")
      nb.anne <- trunc(time)
      dateapdiag <- as.Date(data[, rmap$year] + scale)
      nbj <- difftime(dateapdiag,
                      data[, rmap$year],
                      units = "days")
      coef01 <-
        coef02 <-
        coef03 <-
        coef04 <-
        coef05 <- coef06 <- numeric(length = length(nb.anne))

      id1 <- which(nb.anne != 0 &
                     dateAnniv != data[, rmap$year])
      if (length(id1) != 0) {
        id1AF <- intersect(id1, which(dateAnniv <= as.Date(paste0(
          format(data[, rmap$year], '%Y'), "-12-31"
        ))))


        if (length(id1AF) != 0) {
          coef01[id1AF] <- difftime(dateAnniv[id1AF],
                                    data[id1AF, rmap$year],
                                    units = "days")

          coef02[id1AF] <- difftime(as.Date(paste0(format(data[id1AF, rmap$year], '%Y'), "-12-31")),
                                    dateAnniv[id1AF], units = "days")

          coef03[id1AF] <-
            nbj[id1AF] - coef01[id1AF] - coef02[id1AF]
          id11 <- intersect(id1AF, which(dateDC <= as.Date(paste(
            format(dateAnniv + nbj * nb.anne, '%Y'),
            format(dateAnniv, '%m-%d'),
            sep = '-'
          ))))


          if (length(id11) != 0) {
            res <- data[id11, rmap$year] + nbj[id11] * nb.anne[id11]
            res2 <- try(class(as.Date(res)), TRUE)
            if (inherits(res2, "try-error")) {
              resfin <- paste(format(data[id11, rmap$year] + nbj[id11] * nb.anne[id11], '%Y'),
                              format(data[id11, rmap$year] + 1, '%m-%d'),
                              sep = '-')
            } else{
              resfin <- data[id11, rmap$year] + nbj[id11] * nb.anne[id11]
            }

            coef04[id11] <- difftime(dateDC[id11], as.Date(resfin),
                                     units = "days")
            if (any(coef04[id11] < 0)) {
              coef04[id11][which(coef04[id11] < 0)] <- 0
            }
          }

          id12 <- intersect(id1AF,
                            which(dateDC > (dateAnniv + nbj * nb.anne) &
                              dateDC <= as.Date(paste0(
                                format(data[, rmap$year] + nbj * nb.anne, '%Y'),
                                "-12-31"
                              ))))


          if (length(id12) != 0) {
            coef05[id12] <- difftime(dateDC[id12], as.Date(dateAnniv[id12] + nbj[id12] * nb.anne[id12]), units = "days")
          }
          id13 <- intersect(id1AF, which(dateDC >= as.Date(paste0(
            format(data[, rmap$year] + nbj * (1 + nb.anne), '%Y'), "-01-01"
          ))))


          if (length(id13) != 0) {
            coef06[id13] <- difftime(dateDC[id13],
                                     as.Date(paste0(
                                       format(data[id13, rmap$year] +
                                                nbj[id13] * (1 + nb.anne[id13]), '%Y'),
                                       "-01-01"
                                     )),
                                     units = "days")
          }
        }



        id1FA <- intersect(id1,
                           which(dateAnniv > as.Date(paste0(
                             format(data[, rmap$year], '%Y'), "-12-31"
                           ))))
        if (length(id1FA) != 0) {
          coef01[id1FA] <- difftime(as.Date(paste0(format(data[id1FA, rmap$year], '%Y'), "-12-31")),
                                    data[id1FA, rmap$year], units = "days")

          coef02[id1FA] <- difftime(dateAnniv[id1FA],
                                    as.Date(paste0(
                                      format(data[id1FA, rmap$year] + nbj[id1FA],
                                             '%Y'),
                                      "-01-01"
                                    )),
                                    units = "days")
          coef03[id1FA] <-
            nbj[id1FA] - coef01[id1FA] - coef02[id1FA]
          id14 <- intersect(id1FA, which(dateDC <= as.Date(paste0(
            format(data[, rmap$year] + nbj * nb.anne, '%Y'),
            "-12-31"
          ))))


          if (length(id14) != 0) {
            coef04[id14] <- difftime(dateDC[id14],
                                     c(data[id14, rmap$year] +
                                         nbj[id14] * nb.anne[id14]),
                                     units = "days")
            if(any(coef04[id14] < 0)){
              coef04[id14][which(coef04[id14] < 0)] <- 0
            }

          }
          id15 <- intersect(id1FA,
                            which(dateDC <= as.Date(dateAnniv + nbj * nb.anne) &
                              dateDC >= as.Date(paste0(
                                format(data[, rmap$year] +
                                         nbj * (1 + nb.anne), '%Y'),
                                "-01-01"
                              ))))


          if (length(id15) != 0) {
            coef05[id15] = difftime(dateDC[id15],
                                    as.Date(paste0(
                                      format(data[id15, rmap$year] +
                                               nbj[id15] * (1 + nb.anne[id15]),
                                             '%Y'), "-01-01"
                                    )),
                                    units = "days")
          }

          id16 <- intersect(id1FA,
                            which(dateDC > as.Date(dateAnniv + nbj * nb.anne)))

          if (length(id16) != 0) {
            coef06[id16]  <- difftime(dateDC[id16],
                                      as.Date(dateAnniv[id16] +
                                                nbj[id16] * nb.anne[id16]),
                                      units = "days")
          }
        }
      }

      id2 <- which(nb.anne != 0 &
                     dateAnniv == data[, rmap$year])
      if (length(id2) != 0) {
        coef01[id2] <- difftime(as.Date(paste0(format(data[id2, rmap$year], '%Y'),
                                  "-12-31")),
                   data[id2, rmap$year], units = "days")
        coef02[id2] <- nbj[id2] - coef01[id2]
        #which have died before the yearDC - 12-31
        id21 <- intersect(id2,
                          which(dateDC <= as.Date(paste0(
                            format(data[, rmap$year] + nbj * nb.anne, '%Y'),
                            "-12-31"
                          ))))
        if (length(id21) != 0) {
          coef03[id21] <- difftime(dateDC[id21],
                                   as.Date(data[id21, rmap$year] +
                                             nbj[id21] * nb.anne[id21]),
                                   units = "days")

          coef03[id21][which(coef03[id21] < 0)] <- 0
        }
        id22 <- intersect(id2, which(dateDC >= as.Date(paste0(
          format(data[, rmap$year] + nbj * (1 + nb.anne), '%Y'), "-01-01"
        ))))

        if (length(id22) != 0) {
          coef04[id22] <- difftime(dateDC[id22],
                                   as.Date(paste0(
                                     format(data[id22, rmap$year] +
                                              nbj[id22] * (1 + nb.anne[id22]),
                                            '%Y'),
                                     "-01-01")),
                                   units = "days")
        }

      }

      id3 <- which(nb.anne == 0 & dateAnniv != data[, rmap$year])
      if (length(id3) != 0) {
        id3AF <- intersect(id3,
                           which(dateAnniv <= as.Date(paste0(
                             format(data[, rmap$year], '%Y'), "-12-31"
                           ))))

        if (length(id3AF) != 0) {
          id31 <- intersect(id3AF, which(dateDC <= as.Date(dateAnniv)))

          if (length(id31) != 0) {
            coef01[id31] <- difftime(dateDC[id31],
                                     data[id31, rmap$year],
                                     units = "days")
          }

          id32 <- intersect(id3AF, which(dateDC > as.Date(paste(
            format(dateAnniv, '%Y'),
            format(dateAnniv, '%m-%d'),
            sep = '-'
          )) &
            dateDC <= as.Date(paste0(
              format(data[, rmap$year], '%Y'), "-12-31"
            ))))

          if (length(id32) != 0) {
            coef01[id32] <- difftime(dateAnniv[id32],
                                     data[id32, rmap$year],
                                     units = "days")

            coef02[id32] <- difftime(dateDC[id32],
                                     as.Date(dateAnniv[id32]),
                                     units = "days")
          }

          id33 <- intersect(id3AF,
                            which(dateDC >= as.Date(paste0(
                              format(data[, rmap$year] +
                                       nbj * (1 + nb.anne),
                                     '%Y'),
                              "-01-01"
                            ))))



          if (length(id33) != 0) {
            coef01[id33] <- difftime(dateAnniv[id33],
                                     data[id33, rmap$year],
                                     units = "days")
            coef02[id33] <- difftime(as.Date(paste0(format(
              data[id33, rmap$year], '%Y'
            ), "-12-31")),
            dateAnniv[id33],
            units = "days")

            coef03[id33] <- difftime(dateDC[id33],
                                     as.Date(paste0(
                                       format(data[id33, rmap$year] +
                                                nbj[id33] * (1 + nb.anne[id33]), '%Y'),
                                       "-01-01"
                                     )),
                                     units = "days")

          }
        }

        id3FA <- intersect(id3,
                           which(dateAnniv >= as.Date(paste0(
                             format(data[, rmap$year] + nbj, '%Y'),
                             "-01-01"
                           ))))

        if (length(id3FA) != 0) {
          id34 <- intersect(id3FA,
                            which(dateDC <= as.Date(paste0(
                              format(data[, rmap$year], '%Y'),
                              "-12-31"
                            ))))

          if (length(id34) != 0) {
            coef01[id34] <- difftime(dateDC[id34],
                                     data[id34, rmap$year],
                                     units = "days")
          }


          id35 <- intersect(id3FA,
                            which(dateDC <= as.Date(dateAnniv) &
                              dateDC >= as.Date(paste0(
                                format(data[, rmap$year] +
                                         nbj * (1 + nb.anne), '%Y'),
                                "-01-01"
                              ))))

          if (length(id35) != 0) {
            coef01[id35] <-
              difftime(as.Date(paste0(format(
                data[id35, rmap$year], '%Y'
              ),
              "-12-31")),
              data[id35, rmap$year],
              units = "days")

            coef02[id35] <- difftime(dateDC[id35],
                                     as.Date(paste0(
                                       format(data[id35, rmap$year] +
                                                nbj[id35] * (1 + nb.anne[id35]), '%Y'),
                                       "-01-01"
                                     )),
                                     units = "days")
          }


          id36 <- intersect(id3FA, which(dateDC > as.Date(dateAnniv)))

          if (length(id36) != 0) {
            coef01[id36] <- difftime(as.Date(paste0(format(
              data[id36, rmap$year], '%Y'), "-12-31")),
            data[id36, rmap$year], units = "days")

            coef02[id36] <- difftime(dateAnniv[id36],
                                     as.Date(paste0(
                                       format(data[id36, rmap$year] +
                                                nbj[id36], '%Y'), "-01-01")),
                                     units = "days")

            coef03[id36] <- difftime(dateDC[id36], as.Date(paste(
              format(dateAnniv[id36], '%Y'),
              format(dateAnniv[id36], '%m-%d'),
              sep = '-'
            )),
            units = "days")
          }
        }
      }


      id4 <- which(nb.anne == 0 & dateAnniv == data[, rmap$year])

      if (length(id4) != 0) {
        id41 <- intersect(id4,
                          which(dateDC <= as.Date(paste0(
                            format(data[, rmap$year], '%Y'), "-12-31"
                          ))))


        if (length(id41) != 0) {
          coef01[id41] <- difftime(dateDC[id41],
                                   data[id41, rmap$year],
                                   units = "days")
        }



        id42 <- intersect(id4, which(dateDC >= as.Date(paste0(
          format(data[, rmap$year] + nbj * (1 + nb.anne), '%Y'),
          "-01-01"
        ))))


        if (length(id42) != 0) {
          coef01[id42] <- difftime(as.Date(paste0(format(data[id42, rmap$year], '%Y'), "-12-31")),
                                   data[id42, rmap$year], units = "days")

          coef02[id42] <- difftime(dateDC[id42],
                                   as.Date(paste0(
                                     format(data[id42, rmap$year] +
                                              nbj[id42] * (1 + nb.anne[id42]),
                                            '%Y'), "-01-01"
                                   )),
                                   units = "days")
        }
      }


      coefs <- data.frame(coef01, coef02, coef03, coef04, coef05, coef06)
      dateDiag <- data[, rmap$year]

      ehazardInt <- mapply(
        FUN = cumpop,
        1:nrow(data),
        MoreArgs = list(
          ratetable,
          Indic02,
          nb.anne,
          coefs,
          nb.age.dc,
          dateDiag,
          dateAnniv
        )
      )



      return(list(
        ehazard = ehazard,
        ehazardInt = ehazardInt,
        dateDiag = dateDiag
      ))


    }
  }
}
