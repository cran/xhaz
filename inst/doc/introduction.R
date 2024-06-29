## ----include = FALSE----------------------------------------------------------
options(rmarkdown.html_vignette.check_title = FALSE)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(xhaz)

## -----------------------------------------------------------------------------
library(xhaz)

## -----------------------------------------------------------------------------
data("simuData", package = "xhaz")

head(simuData)
dim(simuData)

interval <- c(0, 0.718, 1.351, 2.143, 3.601, 6)
fit.estv1 <- xhaz(formula = Surv(time_year, status) ~ agec + race,
                  data = simuData,
                  ratetable = survexp.us,
                  interval = interval,
                  rmap = list(age = 'age', sex = 'sex', year = 'date'),
                  baseline = "constant",
                  pophaz = "classic")
                  
fit.estv1

## -----------------------------------------------------------------------------
fit.corrected1 <- xhaz(formula = Surv(time_year, status) ~ agec + race,
                       data = simuData,
                       ratetable = survexp.us,
                       interval = interval,
                       rmap = list(age = 'age', sex = 'sex', year = 'date'),
                       baseline = "constant", pophaz = "corrected",
                       add.rmap = "race")
                       
fit.corrected1

## -----------------------------------------------------------------------------
 # An additionnal cavariate (here race) missing in the life tables is
 # considered by the model with a breakpoint at 75 years

 fit.corrected2 <- xhaz(formula = Surv(time_year, status) ~ agec + race,
                        data = simuData,
                        ratetable = survexp.us,
                        interval = interval,
                        rmap = list(age = 'age', sex = 'sex', year = 'date'),
                        baseline = "constant", pophaz = "corrected",
                        add.rmap = "race",
                        add.rmap.cut = list(breakpoint = TRUE, cut = 75))



 fit.corrected2

## -----------------------------------------------------------------------------
AIC(fit.estv1)
AIC(fit.corrected1)

BIC(fit.estv1)
BIC(fit.corrected1)

## -----------------------------------------------------------------------------
anova(fit.corrected1, fit.corrected2)

## ----fig.width=10, fig.height=10----------------------------------------------

#only add Surv variables (time_year and status) to have them in the new.data. 
#They are not used for the prediction

simuData[1,]


newdata1 <-
  expand.grid(
    race = factor("black",
                   levels = levels(simuData$race)),
    agec  = simuData[1, "agec"], #i.e age 50.5 years
    time_year = 0,
    status = 0
  )



predict_mod1 <- predict(object = fit.estv1,
                        new.data = newdata1,
                        times.pts = c(seq(0, 6, 0.1)),
                        baseline = FALSE)


predict_mod2 <- predict(object = fit.corrected1,
                        new.data = newdata1,
                        times.pts = c(seq(0, 6, 0.1)),
                        baseline = FALSE)



predict_mod3 <- predict(object = fit.corrected2,
                        new.data = newdata1,
                        times.pts = c(seq(0, 6, 0.1)),
                        baseline = FALSE)
old.par <- par(no.readonly = TRUE)
par(mfrow = c(2, 1))


plot(
  predict_mod1,
  what = "survival",
  xlab = "time since diagnosis (year)",
  ylab = "net survival",
  ylim = c(0, 1),
  main = "Estève Model"
)

par(new = TRUE)

plot(
  predict_mod2,
  what = "survival",
  xlab = "",
  ylab = "",
  ylim = c(0, 1),
  lty = 2,
  lwd = 2# main = "Touraine Model"
)

par(new = TRUE)
plot(
  predict_mod3,
  what = "survival",
  xlab = "",
  ylab = "",
  ylim = c(0, 1),
  lty = 3,
  lwd = 3#main = "Mba Model"
)
legend(
  "bottomleft",
  legend = c("Esteve Model",
             "Touraine Model",
             "Mba Model"),
  lty = c(1, 2 , 3),
  lwd = c(1, 2 , 3)
)


plot(
  predict_mod1,
  what = "hazard",
  xlab = "time since diagnosis (year)",
  ylab = "excess hazard",
  ylim = c(0, 0.30),
  lty = 1
)

par(new = TRUE)
plot(
  predict_mod2,
  what = "hazard",
  xlab = "",
  ylab = "",
  ylim = c(0, 0.30),
  lty = 2,
  lwd = 2
)
par(new = TRUE)
plot(
  predict_mod3,
  what = "hazard",
  xlab = "",
  ylab = "",
  ylim = c(0, 0.30),
  lty = 3,
  lwd = 3
)

legend(
  "topright",
  legend = c("Esteve Model",
             "Touraine Model",
             "Mba Model"),
  lty = c(1, 2 , 3),
  lwd = c(1, 2 , 3)
)
par(old.par)


## ----fig.width=10, fig.height=10----------------------------------------------

#only add Surv variables (time_year and status) to have them in the new.data. 
#They are not used for the prediction
#we could be interested to the prediction of net survival and excess hazard of the individual we the same characteristics that this one (age 50.5, race black)

predmar_mod1 <- predict(object = fit.estv1,
                        new.data = simuData,
                        times.pts = c(seq(0, 6, 0.1)),
                        baseline = FALSE)


predmar_mod2 <- predict(object = fit.corrected1,
                        new.data = simuData,
                        times.pts = c(seq(0, 6, 0.1)),
                        baseline = FALSE)



predmar_mod3 <- predict(object = fit.corrected2,
                        new.data = simuData,
                        times.pts = c(seq(0, 6, 0.1)),
                        baseline = FALSE)

par(mfrow = c(2, 1))


plot(
  predmar_mod1,
  what = "survival",
  xlab = "time since diagnosis (year)",
  ylab = "net survival",
  ylim = c(0, 1),
  main = "Estève Model"
)

par(new = TRUE)
plot(
  predmar_mod2,
  what = "survival",
  xlab = "",
  ylab = "",
  ylim = c(0, 1),
  lty = 2,
  lwd = 2# main = "Touraine Model"
)

par(new = TRUE)
plot(
  predmar_mod3,
  what = "survival",
  xlab = "",
  ylab = "",
  ylim = c(0, 1),
  lty = 3,
  lwd = 3#main = "Mba Model"
)
legend(
  "bottomleft",
  legend = c("Esteve Model",
             "Touraine Model",
             "Mba Model"),
  lty = c(1, 2 , 3),
  lwd = c(1, 2 , 3)
)


plot(
  predmar_mod1,
  what = "hazard",
  xlab = "time since diagnosis (year)",
  ylab = "excess hazard",
  ylim = c(0, 0.30),
  lty = 1
)

par(new = TRUE)
plot(
  predmar_mod2,
  what = "hazard",
  xlab = "",
  ylab = "",
  ylim = c(0, 0.30),
  lty = 2,
  lwd = 2
)
par(new = TRUE)
plot(
  predmar_mod3,
  what = "hazard",
  xlab = "",
  ylab = "",
  ylim = c(0, 0.30),
  lty = 3,
  lwd = 3
)

legend(
  "topright",
  legend = c("Esteve Model",
             "Touraine Model",
             "Mba Model"),
  lty = c(1, 2 , 3),
  lwd = c(1, 2 , 3)
)
par(old.par)

