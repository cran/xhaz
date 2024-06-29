

cumpop <- function(i,
                   ratetable,
                   Indic00,
                   nb.anne00,
                   coefs00,
                   nb.age.dc,
                   dateDiag00,
                   dateAnniv00) {
  Indic <- as.data.frame(Indic00)[i,]
  coef01 <- coefs00$coef01[i]
  coef02 <- coefs00$coef02[i]
  coef03 <- coefs00$coef03[i]
  coef04 <- coefs00$coef04[i]
  coef05 <- coefs00$coef05[i]
  coef06 <- coefs00$coef06[i]
  nb.anne <- nb.anne00[i]
  dateDiag <- dateDiag00[i]
  dateAnniv <- dateAnniv00[i]

  if (as.integer(nb.anne) != 0) {
    if (dateAnniv != dateDiag) {
      if (dateAnniv <= as.Date(paste0(format(dateDiag, '%Y'), "-12-31"))) {
        if (coef06 != 0) {
          indic01 <-
            do.call(rbind, lapply(1:(nb.anne + 1), function(x)
              rbind(Indic)))
          indic01$age <- indic01$age + c(0:(nb.anne))
          indic01$year <- indic01$year + c(0:(nb.anne))

          indic02 <-
            do.call(rbind, lapply(1:(nb.anne + 1), function(x)
              rbind(Indic)))
          indic02$age <- indic02$age + c(1:(nb.anne + 1))
          indic02$year <- indic02$year + c(0:nb.anne)

          indic03 <- do.call(rbind, lapply(1:nb.anne, function(x)
            rbind(Indic)))
          indic03$age <- indic03$age + c(1:nb.anne)
          indic03$year <- indic03$year + c(1:nb.anne)

          if (max(indic01$age) < nb.age.dc) {
            indic06 <- indic03[nb.anne,]
            indic06$age <- indic06$age + 1
            indic06$year <- indic06$year + 1
          } else{
            indic01$age[indic01$age >= nb.age.dc] <- nb.age.dc
            indic02$age[indic02$age >= nb.age.dc] <- nb.age.dc
            indic03$age[indic03$age >= nb.age.dc] <- nb.age.dc
            indic06 <- indic03[nb.anne,]
            indic06$age <- indic06$age + 1
            indic06$year <- indic06$year + 1
          }


          if (any(indic01$age>
                  max(order(attributes(ratetable)$dimnames[[1]])))) {
            indic01$age[
              which(indic01$age>
                      max(order(attributes(ratetable)$dimnames[[1]])))
            ] <- max(order(attributes(ratetable)$dimnames[[1]]))
          }

          if (any(indic02$age>
                  max(order(attributes(ratetable)$dimnames[[1]])))) {
            indic02$age[
              which(indic02$age>
                      max(order(attributes(ratetable)$dimnames[[1]])))
            ] <- max(order(attributes(ratetable)$dimnames[[1]]))
          }

          if (any(indic03$age>
                  max(order(attributes(ratetable)$dimnames[[1]])))) {
            indic03$age[
              which(indic03$age>
                      max(order(attributes(ratetable)$dimnames[[1]])))
            ] <- max(order(attributes(ratetable)$dimnames[[1]]))
          }

          if (any(indic06$age>
                  max(order(attributes(ratetable)$dimnames[[1]])))) {
            indic06$age[
              which(indic06$age>
                      max(order(attributes(ratetable)$dimnames[[1]])))
            ] <- max(order(attributes(ratetable)$dimnames[[1]]))
          }


          results <- sum(coef01 * ratetable[as.matrix(indic01)]) +
            sum(coef02 * ratetable[as.matrix(indic02)]) +
            sum(coef03 * ratetable[as.matrix(indic03)]) +
            coef06 * ratetable[as.matrix(indic06)]
        }
        else{
          if (coef05 != 0) {
            indic01 <- do.call(rbind,
                               lapply(1:(nb.anne + 1),
                                      function(x)
                                        rbind(Indic)))
            indic01$age <- indic01$age + c(0:(nb.anne))
            indic01$year <- indic01$year + c(0:(nb.anne))

            indic02 <- do.call(rbind,
                               lapply(1:nb.anne,
                                      function(x)
                                        rbind(Indic)))
            indic02$age <- indic02$age + c(1:nb.anne)
            indic02$year <- indic02$year + c(0:(nb.anne - 1))

            indic03 <- do.call(rbind,
                               lapply(1:nb.anne,
                                      function(x)
                                        rbind(Indic)))

            indic03$age <- indic03$age + c(1:nb.anne)
            indic03$year <- indic03$year + c(1:nb.anne)

            if (max(indic01$age) < nb.age.dc) {
              indic05 <- indic02[nb.anne,]
              indic05$age <- indic05$age + 1
              indic05$year <- indic05$year + 1
            } else{
              indic01$age[indic01$age >= nb.age.dc] <- nb.age.dc
              indic02$age[indic02$age >= nb.age.dc] <- nb.age.dc
              indic03$age[indic03$age >= nb.age.dc] <- nb.age.dc
              indic05 <- indic02[nb.anne,]
              indic05$age <- indic05$age + 1
              indic05$year <- indic05$year + 1
            }

            if (any(indic01$age>
                    max(order(attributes(ratetable)$dimnames[[1]])))) {
              indic01$age[
                which(indic01$age>
                        max(order(attributes(ratetable)$dimnames[[1]])))
              ] <- max(order(attributes(ratetable)$dimnames[[1]]))
            }

            if (any(indic02$age>
                    max(order(attributes(ratetable)$dimnames[[1]])))) {
              indic02$age[
                which(indic02$age>
                        max(order(attributes(ratetable)$dimnames[[1]])))
              ] <- max(order(attributes(ratetable)$dimnames[[1]]))
            }

            if (any(indic03$age>
                    max(order(attributes(ratetable)$dimnames[[1]])))) {
              indic03$age[
                which(indic03$age>
                        max(order(attributes(ratetable)$dimnames[[1]])))
              ] <- max(order(attributes(ratetable)$dimnames[[1]]))
            }

            if (any(indic05$age>
                    max(order(attributes(ratetable)$dimnames[[1]])))) {
              indic05$age[
                which(indic05$age>
                        max(order(attributes(ratetable)$dimnames[[1]])))
              ] <- max(order(attributes(ratetable)$dimnames[[1]]))
            }
            results <- sum(coef01 * ratetable[as.matrix(indic01)]) +
              sum(coef02 * ratetable[as.matrix(indic02)]) +
              sum(coef03 * ratetable[as.matrix(indic03)]) +
              coef05 * ratetable[as.matrix(indic05)]
          }
          else{
            if (coef04 != 0) {
              indic01 <- do.call(rbind,
                                 lapply(1:nb.anne,
                                        function(x)
                                          rbind(Indic)))
              indic01$age <- indic01$age + c(0:(nb.anne - 1))
              indic01$year <- indic01$year + c(0:(nb.anne - 1))

              indic02 <-
                do.call(rbind, lapply(1:nb.anne, function(x)
                  rbind(Indic)))
              indic02$age <- indic02$age + c(1:nb.anne)
              indic02$year <- indic02$year + c(0:(nb.anne - 1))
              indic03 <- do.call(rbind,
                                 lapply(1:nb.anne, function(x)rbind(Indic)))
              indic03$age <- indic03$age + c(1:nb.anne)
              indic03$year <- indic03$year + c(1:nb.anne)

              if (max(indic01$age) < nb.age.dc) {
                indic04 <- indic01[nb.anne,]
                indic04$age <- indic04$age + 1
                indic04$year <- indic04$year + 1
              } else{
                indic01$age[indic01$age >= nb.age.dc] <- nb.age.dc
                indic02$age[indic02$age >= nb.age.dc] <- nb.age.dc
                indic03$age[indic03$age >= nb.age.dc] <- nb.age.dc
                indic04 <- indic01[nb.anne,]
                indic04$age <- indic04$age + 1
                indic04$year <- indic04$year + 1
              }

              if (any(indic01$age>
                      max(order(attributes(ratetable)$dimnames[[1]])))) {
                indic01$age[
                  which(indic01$age>
                          max(order(attributes(ratetable)$dimnames[[1]])))
                ] <- max(order(attributes(ratetable)$dimnames[[1]]))
              }

              if (any(indic02$age>
                      max(order(attributes(ratetable)$dimnames[[1]])))) {
                indic02$age[
                  which(indic02$age>
                          max(order(attributes(ratetable)$dimnames[[1]])))
                ] <- max(order(attributes(ratetable)$dimnames[[1]]))
              }

              if (any(indic03$age>
                      max(order(attributes(ratetable)$dimnames[[1]])))) {
                indic03$age[
                  which(indic03$age>
                          max(order(attributes(ratetable)$dimnames[[1]])))
                ] <- max(order(attributes(ratetable)$dimnames[[1]]))
              }

              if (any(indic04$age>
                      max(order(attributes(ratetable)$dimnames[[1]])))) {
                indic04$age[
                  which(indic04$age>
                          max(order(attributes(ratetable)$dimnames[[1]])))
                ] <- max(order(attributes(ratetable)$dimnames[[1]]))
              }
              results <- sum(coef01 * ratetable[as.matrix(indic01)]) +
                sum(coef02 * ratetable[as.matrix(indic02)]) +
                sum(coef03 * ratetable[as.matrix(indic03)]) +
                coef04 * ratetable[as.matrix(indic04)]
            }
            else{
              indic01 <-
                do.call(rbind, lapply(1:nb.anne, function(x)
                  rbind(Indic)))
              indic01$age <- indic01$age + c(0:(nb.anne - 1))
              indic01$year <- indic01$year + c(0:(nb.anne - 1))

              indic02 <-
                do.call(rbind, lapply(1:nb.anne, function(x)
                  rbind(Indic)))
              indic02$age <- indic02$age + c(1:nb.anne)
              indic02$year <- indic02$year + c(0:(nb.anne - 1))

              indic03 <-
                do.call(rbind, lapply(1:nb.anne, function(x)
                  rbind(Indic)))
              indic03$age <- indic03$age + c(1:nb.anne)
              indic03$year <- indic03$year + c(1:nb.anne)


              if (any(indic01$age>
                      max(order(attributes(ratetable)$dimnames[[1]])))) {
                indic01$age[
                  which(indic01$age>
                          max(order(attributes(ratetable)$dimnames[[1]])))
                ] <- max(order(attributes(ratetable)$dimnames[[1]]))
              }

              if (any(indic02$age>
                      max(order(attributes(ratetable)$dimnames[[1]])))) {
                indic02$age[
                  which(indic02$age>
                          max(order(attributes(ratetable)$dimnames[[1]])))
                ] <- max(order(attributes(ratetable)$dimnames[[1]]))
              }

              if (any(indic03$age>
                      max(order(attributes(ratetable)$dimnames[[1]])))) {
                indic03$age[
                  which(indic03$age>
                          max(order(attributes(ratetable)$dimnames[[1]])))
                ] <- max(order(attributes(ratetable)$dimnames[[1]]))
              }


              results <-
                sum(coef01 * ratetable[as.matrix(indic01)]) +
                sum(coef02 * ratetable[as.matrix(indic02)]) +
                sum(coef03 * ratetable[as.matrix(indic03)])
            }
          }
        }
      }
      else{
        if (coef06 != 0) {
          indic01 <-
            do.call(rbind, lapply(1:(nb.anne + 1), function(x)
              rbind(Indic)))
          indic01$age <- indic01$age + c(0:(nb.anne))
          indic01$year <- indic01$year + c(0:(nb.anne))

          indic02 <- do.call(rbind, lapply(1:(nb.anne + 1), function(x)
            rbind(Indic)))
          indic02$age <- indic02$age + c(0:nb.anne)
          indic02$year <- indic02$year + c(1:(nb.anne + 1))

          indic03 <- do.call(rbind, lapply(1:nb.anne, function(x)
            rbind(Indic)))
          indic03$age <- indic03$age + c(1:nb.anne)
          indic03$year <- indic03$year + c(1:nb.anne)

          if (max(indic01$age) < nb.age.dc) {
            indic06 <- indic03[nb.anne,]
            indic06$age <- indic06$age + 1
            indic06$year <- indic06$year + 1
          } else{
            indic01$age[indic01$age >= nb.age.dc] <- nb.age.dc
            indic02$age[indic02$age >= nb.age.dc] <- nb.age.dc
            indic03$age[indic03$age >= nb.age.dc] <- nb.age.dc
            indic06 <- indic03[nb.anne,]
            indic06$age <- indic06$age + 1
            indic06$year <- indic06$year + 1
          }

          if (any(indic01$age>
                  max(order(attributes(ratetable)$dimnames[[1]])))) {
            indic01$age[
              which(indic01$age>
                      max(order(attributes(ratetable)$dimnames[[1]])))
            ] <- max(order(attributes(ratetable)$dimnames[[1]]))
          }

          if (any(indic02$age>
                  max(order(attributes(ratetable)$dimnames[[1]])))) {
            indic02$age[
              which(indic02$age>
                      max(order(attributes(ratetable)$dimnames[[1]])))
            ] <- max(order(attributes(ratetable)$dimnames[[1]]))
          }

          if (any(indic03$age>
                  max(order(attributes(ratetable)$dimnames[[1]])))) {
            indic03$age[
              which(indic03$age>
                      max(order(attributes(ratetable)$dimnames[[1]])))
            ] <- max(order(attributes(ratetable)$dimnames[[1]]))
          }

          if (any(indic06$age>
                  max(order(attributes(ratetable)$dimnames[[1]])))) {
            indic06$age[
              which(indic06$age>
                      max(order(attributes(ratetable)$dimnames[[1]])))
            ] <- max(order(attributes(ratetable)$dimnames[[1]]))
          }
          results <- sum(coef01 * ratetable[as.matrix(indic01)]) +
            sum(coef02 * ratetable[as.matrix(indic02)]) +
            sum(coef03 * ratetable[as.matrix(indic03)]) +
            coef06 * ratetable[as.matrix(indic06)]
          }
        else{
          if (coef05 != 0) {
            indic01 <- do.call(rbind, lapply(1:(nb.anne + 1), function(x)
              rbind(Indic)))
            indic01$age <- indic01$age + c(0:(nb.anne))
            indic01$year <- indic01$year + c(0:(nb.anne))

            indic02 <- do.call(rbind, lapply(1:nb.anne, function(x)
              rbind(Indic)))
            indic02$age <- indic02$age + c(0:(nb.anne - 1))
            indic02$year <- indic02$year + c(1:nb.anne)

            indic03 <- do.call(rbind, lapply(1:nb.anne, function(x)
              rbind(Indic)))
            indic03$age <- indic03$age + c(1:nb.anne)
            indic03$year <- indic03$year + c(1:nb.anne)

            if (max(indic01$age) < nb.age.dc) {
              indic05 <- indic02[nb.anne,]
              indic05$age <- indic05$age + 1
              indic05$year <- indic05$year + 1
            } else{
              indic01$age[indic01$age >= nb.age.dc] <- nb.age.dc
              indic02$age[indic02$age >= nb.age.dc] <- nb.age.dc
              indic03$age[indic03$age >= nb.age.dc] <- nb.age.dc
              indic05 <- indic02[nb.anne,]
              indic05$age <- indic05$age + 1
              indic05$year <- indic05$year + 1
            }


            if (any(indic01$age>
                    max(order(attributes(ratetable)$dimnames[[1]])))) {
              indic01$age[
                which(indic01$age>
                        max(order(attributes(ratetable)$dimnames[[1]])))
              ] <- max(order(attributes(ratetable)$dimnames[[1]]))
            }

            if (any(indic02$age>
                    max(order(attributes(ratetable)$dimnames[[1]])))) {
              indic02$age[
                which(indic02$age>
                        max(order(attributes(ratetable)$dimnames[[1]])))
              ] <- max(order(attributes(ratetable)$dimnames[[1]]))
            }

            if (any(indic03$age>
                    max(order(attributes(ratetable)$dimnames[[1]])))) {
              indic03$age[
                which(indic03$age>
                        max(order(attributes(ratetable)$dimnames[[1]])))
              ] <- max(order(attributes(ratetable)$dimnames[[1]]))
            }

            if (any(indic05$age>
                    max(order(attributes(ratetable)$dimnames[[1]])))) {
              indic05$age[
                which(indic05$age>
                        max(order(attributes(ratetable)$dimnames[[1]])))
              ] <- max(order(attributes(ratetable)$dimnames[[1]]))
            }
            results <-
              sum(coef01 * ratetable[as.matrix(indic01)]) +
              sum(coef02 * ratetable[as.matrix(indic02)]) +
              sum(coef03 * ratetable[as.matrix(indic03)]) +
              coef05 * ratetable[as.matrix(indic05)]
          } else{

            if (coef04 != 0) {
              indic01 <- do.call(rbind, lapply(1:nb.anne, function(x)
                rbind(Indic)))
              indic01$age <- indic01$age + c(0:(nb.anne - 1))
              indic01$year <- indic01$year + c(0:(nb.anne - 1))

              indic02 <- do.call(rbind, lapply(1:nb.anne, function(x)
                rbind(Indic)))
              indic02$age <- indic02$age + c(0:(nb.anne - 1))
              indic02$year <- indic02$year + c(1:nb.anne)

              indic03 <- do.call(rbind, lapply(1:nb.anne, function(x)
                rbind(Indic)))
              indic03$age <- indic03$age + c(1:nb.anne)
              indic03$year <- indic03$year + c(1:nb.anne)

              if (max(indic01$age) < nb.age.dc) {
                indic04 <- indic01[nb.anne,]
                indic04$age <- indic04$age + 1
                indic04$year <- indic04$year + 1
              } else{
                indic01$age[indic01$age >= nb.age.dc] <- nb.age.dc
                indic02$age[indic02$age >= nb.age.dc] <- nb.age.dc
                indic03$age[indic03$age >= nb.age.dc] <- nb.age.dc
                indic04 <- indic01[nb.anne,]
                indic04$age <- indic04$age + 1
                indic04$year <- indic04$year + 1
              }



              if (any(indic01$age>
                      max(order(attributes(ratetable)$dimnames[[1]])))) {
                indic01$age[
                  which(indic01$age>
                          max(order(attributes(ratetable)$dimnames[[1]])))
                ] <- max(order(attributes(ratetable)$dimnames[[1]]))
              }

              if (any(indic02$age>
                      max(order(attributes(ratetable)$dimnames[[1]])))) {
                indic02$age[
                  which(indic02$age>
                          max(order(attributes(ratetable)$dimnames[[1]])))
                ] <- max(order(attributes(ratetable)$dimnames[[1]]))
              }

              if (any(indic03$age>
                      max(order(attributes(ratetable)$dimnames[[1]])))) {
                indic03$age[
                  which(indic03$age>
                          max(order(attributes(ratetable)$dimnames[[1]])))
                ] <- max(order(attributes(ratetable)$dimnames[[1]]))
              }

              if (any(indic04$age>
                      max(order(attributes(ratetable)$dimnames[[1]])))) {
                indic04$age[
                  which(indic04$age>
                          max(order(attributes(ratetable)$dimnames[[1]])))
                ] <- max(order(attributes(ratetable)$dimnames[[1]]))
              }

              results <- sum(coef01 * ratetable[as.matrix(indic01)]) +
                sum(coef02 * ratetable[as.matrix(indic02)]) +
                sum(coef03 * ratetable[as.matrix(indic03)]) +
                coef04 * ratetable[as.matrix(indic04)]
            } else{
              indic01 <- do.call(rbind,
                                 lapply(1:nb.anne,
                                        function(x)
                                          rbind(Indic)))

              indic01$age <- indic01$age + c(0:(nb.anne - 1))
              indic01$year <- indic01$year + c(0:(nb.anne - 1))

              indic02 <- do.call(rbind, lapply(1:nb.anne, function(x)
                rbind(Indic)))
              indic02$age <- indic02$age + c(0:(nb.anne - 1))
              indic02$year <- indic02$year + c(1:nb.anne)

              indic03 <- do.call(rbind, lapply(1:nb.anne, function(x)
                rbind(Indic)))
              indic03$age <- indic03$age + c(1:nb.anne)
              indic03$year <- indic03$year + c(1:nb.anne)

              if (any(indic01$age>
                      max(order(attributes(ratetable)$dimnames[[1]])))) {
                indic01$age[
                  which(indic01$age>
                          max(order(attributes(ratetable)$dimnames[[1]])))
                ] <- max(order(attributes(ratetable)$dimnames[[1]]))
              }

              if (any(indic02$age>
                      max(order(attributes(ratetable)$dimnames[[1]])))) {
                indic02$age[
                  which(indic02$age>
                          max(order(attributes(ratetable)$dimnames[[1]])))
                ] <- max(order(attributes(ratetable)$dimnames[[1]]))
              }

              if (any(indic03$age>
                      max(order(attributes(ratetable)$dimnames[[1]])))) {
                indic03$age[
                  which(indic03$age>
                          max(order(attributes(ratetable)$dimnames[[1]])))
                ] <- max(order(attributes(ratetable)$dimnames[[1]]))
              }

              results <- sum(coef01 * ratetable[as.matrix(indic01)]) +
                sum(coef02 * ratetable[as.matrix(indic02)]) +
                sum(coef03 * ratetable[as.matrix(indic03)])
            }
          }
        }
      }
    } else{
      if (coef04 == 0) {
        indic01 <- do.call(rbind, lapply(1:nb.anne, function(x)
          rbind(Indic)))
        indic01$age <- indic01$age + c(0:(nb.anne - 1))
        indic01$year <- indic01$year + c(0:(nb.anne - 1))

        indic02 <- do.call(rbind, lapply(1:nb.anne, function(x) rbind(Indic)))
        indic02$age <- indic02$age + c(0:(nb.anne - 1))
        indic02$year <- indic02$year + c(1:(nb.anne))

        if (max(indic01$age) < nb.age.dc) {
          indic03 <- indic01[nb.anne,]
          indic03$age <- indic03$age + 1
          indic03$year <- indic03$year + 1
        } else{
          indic01$age[indic01$age >= nb.age.dc] <- nb.age.dc
          indic02$age[indic02$age >= nb.age.dc] <- nb.age.dc
          indic03 <- indic01[nb.anne,]
          indic03$age <- indic03$age + 1
          indic03$year <- indic03$year + 1
        }


        if (any(indic01$age>
                max(order(attributes(ratetable)$dimnames[[1]])))) {
          indic01$age[
            which(indic01$age>
                    max(order(attributes(ratetable)$dimnames[[1]])))
          ] <- max(order(attributes(ratetable)$dimnames[[1]]))
        }

        if (any(indic02$age>
                max(order(attributes(ratetable)$dimnames[[1]])))) {
          indic02$age[
            which(indic02$age>
                    max(order(attributes(ratetable)$dimnames[[1]])))
          ] <- max(order(attributes(ratetable)$dimnames[[1]]))
        }

        if (any(indic03$age>
                max(order(attributes(ratetable)$dimnames[[1]])))) {
          indic03$age[
            which(indic03$age>
                    max(order(attributes(ratetable)$dimnames[[1]])))
          ] <- max(order(attributes(ratetable)$dimnames[[1]]))
        }



        results <-
          sum(coef01 * ratetable[as.matrix(indic01)]) +
          sum(coef02 * ratetable[as.matrix(indic02)]) +
          coef03 * ratetable[as.matrix(indic03)]
      } else{
        indic01 <- do.call(rbind, lapply(1:(nb.anne + 1), function(x)
            rbind(Indic)))
        indic01$age <- indic01$age + c(0:(nb.anne))
        indic01$year <- indic01$year + c(0:(nb.anne))

        indic02 <- do.call(rbind, lapply(1:nb.anne, function(x)
            rbind(Indic)))
        indic02$age <- indic02$age + c(0:(nb.anne - 1))
        indic02$year <- indic02$year + c(1:(nb.anne))

        if (max(indic01$age) < nb.age.dc) {
          indic04 <- indic02[nb.anne,]
          indic04$age <- indic04$age + 1
          indic04$year <- indic04$year + 1
        } else{
          indic01$age[indic01$age >= nb.age.dc] <- nb.age.dc
          indic02$age[indic02$age >= nb.age.dc] <- nb.age.dc
          indic04 <- indic02[nb.anne,]
          indic04$age <- indic04$age + 1
          indic04$year <- indic04$year + 1
        }

        if (any(indic01$age>
                max(order(attributes(ratetable)$dimnames[[1]])))) {
          indic01$age[
            which(indic01$age>
                    max(order(attributes(ratetable)$dimnames[[1]])))
          ] <- max(order(attributes(ratetable)$dimnames[[1]]))
        }

        if (any(indic02$age>
                max(order(attributes(ratetable)$dimnames[[1]])))) {
          indic02$age[
            which(indic02$age>
                    max(order(attributes(ratetable)$dimnames[[1]])))
          ] <- max(order(attributes(ratetable)$dimnames[[1]]))
        }

        if (any(indic04$age>
                max(order(attributes(ratetable)$dimnames[[1]])))) {
          indic04$age[
            which(indic04$age>
                    max(order(attributes(ratetable)$dimnames[[1]])))
          ] <- max(order(attributes(ratetable)$dimnames[[1]]))
        }

        results <- sum(coef01 * ratetable[as.matrix(indic01)]) +
          sum(coef02 * ratetable[as.matrix(indic02)]) +
          coef04 * ratetable[as.matrix(indic04)]
      }
    }
  } else{
    if (dateAnniv != dateDiag) {
      if (dateAnniv <= as.Date(paste0(format(dateDiag, '%Y'), "-12-31"))) {
        if (coef03 != 0) {
          indic01 <- Indic

          indic02 <- Indic
          indic02$age <- indic02$age + 1

          indic03 <- Indic
          indic03$age <- indic03$age + 1
          indic03$year <- indic03$year + 1


          if (any(indic01$age>
                  max(order(attributes(ratetable)$dimnames[[1]])))) {
            indic01$age[
              which(indic01$age>
                      max(order(attributes(ratetable)$dimnames[[1]])))
            ] <- max(order(attributes(ratetable)$dimnames[[1]]))
          }

          if (any(indic02$age>
                  max(order(attributes(ratetable)$dimnames[[1]])))) {
            indic02$age[
              which(indic02$age>
                      max(order(attributes(ratetable)$dimnames[[1]])))
            ] <- max(order(attributes(ratetable)$dimnames[[1]]))
          }

          if (any(indic03$age>
                  max(order(attributes(ratetable)$dimnames[[1]])))) {
            indic03$age[
              which(indic03$age>
                      max(order(attributes(ratetable)$dimnames[[1]])))
            ] <- max(order(attributes(ratetable)$dimnames[[1]]))
          }



          results <- coef01 * ratetable[as.matrix(indic01)] +
            coef02 * ratetable[as.matrix(indic02)] +
            coef03 * ratetable[as.matrix(indic03)]
          }
        else{
          if (coef02 != 0) {
            indic01 <- Indic

            indic02 <- Indic
            indic02$age <- indic02$age + 1

            if (any(indic01$age>
                    max(order(attributes(ratetable)$dimnames[[1]])))) {
              indic01$age[
                which(indic01$age>
                        max(order(attributes(ratetable)$dimnames[[1]])))
              ] <- max(order(attributes(ratetable)$dimnames[[1]]))
            }

            if (any(indic02$age>
                    max(order(attributes(ratetable)$dimnames[[1]])))) {
              indic02$age[
                which(indic02$age>
                        max(order(attributes(ratetable)$dimnames[[1]])))
              ] <- max(order(attributes(ratetable)$dimnames[[1]]))
            }

            results <- coef01 * ratetable[as.matrix(indic01)] +
              coef02 * ratetable[as.matrix(indic02)]
          } else{
            indic01 <- Indic


            if (any(indic01$age>
                    max(order(attributes(ratetable)$dimnames[[1]])))) {
              indic01$age[
                which(indic01$age>
                        max(order(attributes(ratetable)$dimnames[[1]])))
              ] <- max(order(attributes(ratetable)$dimnames[[1]]))
            }

             results <- coef01 * ratetable[as.matrix(indic01)]
          }
        }
      } else{
        if (coef03 != 0) {
          indic01 <- Indic

          indic02 <- Indic
          indic02$year <- indic02$year + 1

          indic03 <- Indic
          indic03$age <- indic03$age + 1
          indic03$year <- indic03$year + 1

          if (any(indic01$age>
                  max(order(attributes(ratetable)$dimnames[[1]])))) {
            indic01$age[
              which(indic01$age>
                      max(order(attributes(ratetable)$dimnames[[1]])))
            ] <- max(order(attributes(ratetable)$dimnames[[1]]))
          }

          if (any(indic02$age>
                  max(order(attributes(ratetable)$dimnames[[1]])))) {
            indic02$age[
              which(indic02$age>
                      max(order(attributes(ratetable)$dimnames[[1]])))
            ] <- max(order(attributes(ratetable)$dimnames[[1]]))
          }

          if (any(indic03$age>
                  max(order(attributes(ratetable)$dimnames[[1]])))) {
            indic03$age[
              which(indic03$age>
                      max(order(attributes(ratetable)$dimnames[[1]])))
            ] <- max(order(attributes(ratetable)$dimnames[[1]]))
          }



          results <- coef01 * ratetable[as.matrix(indic01)] +
            coef02 * ratetable[as.matrix(indic02)] +
            coef03 * ratetable[as.matrix(indic03)]

        }
        else{
          if (coef02 != 0) {
            indic01 <- Indic

            indic02 <- Indic
            indic02$year <- indic02$year + 1

            if (any(indic01$age>
                    max(order(attributes(ratetable)$dimnames[[1]])))) {
              indic01$age[
                which(indic01$age>
                        max(order(attributes(ratetable)$dimnames[[1]])))
              ] <- max(order(attributes(ratetable)$dimnames[[1]]))
            }

            if (any(indic02$age>
                    max(order(attributes(ratetable)$dimnames[[1]])))) {
              indic02$age[
                which(indic02$age>
                        max(order(attributes(ratetable)$dimnames[[1]])))
              ] <- max(order(attributes(ratetable)$dimnames[[1]]))
            }

            results <- coef01 * ratetable[as.matrix(indic01)] +
              coef02 * ratetable[as.matrix(indic02)]
          }
          else{
            indic01 <- Indic

            if (any(indic01$age>
                    max(order(attributes(ratetable)$dimnames[[1]])))) {
              indic01$age[
                which(indic01$age>
                        max(order(attributes(ratetable)$dimnames[[1]])))
              ] <- max(order(attributes(ratetable)$dimnames[[1]]))
            }

            results <- coef01 * ratetable[as.matrix(indic01)]
          }
        }
      }
    } else{
      if (coef02 != 0) {
        indic01 <- Indic

        indic02 <- Indic
        indic02$year <- indic02$year + 1

        if (any(indic01$age>
                max(order(attributes(ratetable)$dimnames[[1]])))) {
          indic01$age[
            which(indic01$age>
                    max(order(attributes(ratetable)$dimnames[[1]])))
          ] <- max(order(attributes(ratetable)$dimnames[[1]]))
        }

        if (any(indic02$age>
                max(order(attributes(ratetable)$dimnames[[1]])))) {
          indic02$age[
            which(indic02$age>
                    max(order(attributes(ratetable)$dimnames[[1]])))
          ] <- max(order(attributes(ratetable)$dimnames[[1]]))
        }

        results <- coef01 * ratetable[as.matrix(indic01)] +
          coef02 * ratetable[as.matrix(indic02)]
      } else{
        indic01 <- Indic

        if (any(indic01$age>
                max(order(attributes(ratetable)$dimnames[[1]])))) {
          indic01$age[
            which(indic01$age>
                    max(order(attributes(ratetable)$dimnames[[1]])))
          ] <- max(order(attributes(ratetable)$dimnames[[1]]))
        }

      results <- coef01 * ratetable[as.matrix(indic01)]
      }
    }
  }

  return(as.numeric(results))
  }
