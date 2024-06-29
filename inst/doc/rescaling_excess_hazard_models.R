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
data("breast", package = "xhaz")

head(breast)
dim(breast)

# The life table to be used is survexp.us. Note that SEX is coded 2 instead of female in survexp.us.

breast$sexe <- "female"
fit.haz <- exphaz(
               formula = Surv(temps, statut) ~ 1,
               data = breast, ratetable = survexp.us, 
               only_ehazard = FALSE,
               rmap = list(age = 'age', sex = 'sexe',
                           year = 'date'))



breast$expected <- fit.haz$ehazard
breast$expectedCum <- fit.haz$ehazardInt


qknots <- quantile(breast[breast$statut==1,]$temps, probs=c(1:2/3))
mod.bs <- mexhazLT(
  formula = Surv(temps, statut) ~ agecr + armt,
  data = breast, ratetable = survexp.us, degree = 3,
  knots = qknots, expected = "expected", 
  expectedCum = "expectedCum",
  base = "exp.bs", pophaz = "classic")
                  
mod.bs

## -----------------------------------------------------------------------------
mod.bs2 <- mexhazLT(
  formula = Surv(temps, statut) ~ agecr + armt,
  data = breast, ratetable = survexp.us, degree = 3,
  knots = qknots, expected = "expected",
  expectedCum = "expectedCum",
  base = "exp.bs", pophaz = "rescaled")
                  
mod.bs2

## -----------------------------------------------------------------------------

mod.bs3 <- mexhazLT(
  formula = Surv(temps, statut) ~ agecr + armt,
  data = breast, ratetable = survexp.us, degree = 3,
  knots = qknots, expected = "expected",
  expectedCum = "expectedCum",
  base = "exp.bs", pophaz = "classic", random = "hosp")
                  
mod.bs3

## -----------------------------------------------------------------------------

mod.bs4 <- mexhazLT(
  formula = Surv(temps, statut) ~ agecr + armt,
  data = breast, ratetable = survexp.us, degree = 3,
  knots = qknots, expected = "expected", 
  expectedCum = "expectedCum",
  base = "exp.bs", pophaz = "rescaled", random = "hosp")
                  
mod.bs4

## -----------------------------------------------------------------------------
compared_models <- list(mod.bs,mod.bs2, mod.bs3, mod.bs4)
names(compared_models) <- c("mod.bs","mod.bs2", "mod.bs3", "mod.bs4")

sapply(compared_models, function(i) {
  AIC(i)
})


## -----------------------------------------------------------------------------
anova(mod.bs,mod.bs2)
anova(mod.bs3,mod.bs4)

## ----fig.width=10, fig.height=10----------------------------------------------

predict_mod <- predict(mod.bs,
                        time.pts=seq(0.1,10,by=0.1),
                        data.val=data.frame(agecr = breast[1,]$agecr,
                                            armt = breast[1,]$armt))


predict_mod2 <- predict(mod.bs2,
                        time.pts=seq(0.1,10,by=0.1),
                        data.val=data.frame(agecr = breast[1,]$agecr,
                                            armt = breast[1,]$armt))

predict_mod3 <- predict(mod.bs3,
                        time.pts=seq(0.1,10,by=0.1),
                        data.val=data.frame(agecr = breast[1,]$agecr,
                                            armt = breast[1,]$armt, 
                                            hosp = breast[1,]$hosp))

predict_mod4 <- predict(mod.bs4,
                        time.pts=seq(0.1,10,by=0.1),
                        data.val=data.frame(agecr = breast[1,]$agecr,
                                            armt = breast[1,]$armt,
                                            hosp = breast[1,]$hosp))

old.par <- par(no.readonly = TRUE)

par(mfrow = c(1, 2))

plot(predict_mod$results$time.pts, predict_mod$results$hazard,
     type = "l", lwd = 2, xlab = 'Time (years)',
     ylab = "excess hazard")
lines(predict_mod2$results$time.pts, predict_mod2$results$hazard,
      type = "l", lwd = 2, col = "blue", lty = 2)
lines(predict_mod3$results$time.pts, predict_mod3$results$hazard, 
      type = "l", lwd = 2, col = "red")
lines(predict_mod4$results$time.pts, predict_mod4$results$hazard,
      type = "l", lwd = 2, col = "green", lty = 2)


plot(predict_mod$results$time.pts, predict_mod$results$surv,
     type = "l", lwd = 2, xlab = 'Time (years)',
     ylab = "Net survival", ylim = c(0,1))
lines(predict_mod2$results$time.pts, predict_mod2$results$surv,
      type = "l", lwd = 2, col = "blue", lty = 2)
lines(predict_mod3$results$time.pts, predict_mod3$results$surv,
      type = "l", lwd = 2, col = "red")
lines(predict_mod4$results$time.pts, predict_mod4$results$surv,
      type = "l", lwd = 2, col = "green",lty = 2)

legend(
  "topright",
  legend = c("mod.bs",
             "mod.bs2",
             "mod.bs3",
             "mod.bs4"),
  lty = c(1, 2 , 1, 2),
  lwd = c(2, 2 , 2, 2),
  col = c("black", "blue", "red", "green")
)
par(old.par)


