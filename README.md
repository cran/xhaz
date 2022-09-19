
# xhaz: Excess Hazard Modelling Considering Inappropriate Mortality Rates

An R package to fit excess hazard models, with or without proportional
population hazards assumption. The baseline excess hazard could be a
piecewise constant function or a B-splines. When B-splines is choosen
for the baseline excess hazard, the user can specify some covariates
which have a time-dependent effect (using “bsplines”) on the baseline
excess hazard. The user can also specify if the framework corresponds to
the classical excess hazard modeling, i.e. assuming that the expected
mortality of studied individuals is appropriate. He can also consider
two other framework: first, the expected mortality available in the life
table is not accurate and requires taking into account an additional
variable in the life table with a proportional [Touraine et
al. (2020)](https://pubmed.ncbi.nlm.nih.gov/30674229/) or
non-proportional [(Mba et
al. (2020)](https://pubmed.ncbi.nlm.nih.gov/33121436/) effect; second,
there is a non-comparability source of bias in terms of expected
mortality of selected individuals in a clinical trials [Goungounga et
al. (2019)](https://pubmed.ncbi.nlm.nih.gov/31096911/).

(*Here’s the abstract from Touraine et al. paper:* Relative survival
methods used to estimate the excess mortality of cancer patients rely on
the background (or expected) mortality derived from general population
life tables. These methods are based on splitting the observed mortality
into the excess mortality and the background mortality. By assuming a
regression model for the excess mortality, usually a Cox-type model, one
may investigate the effects of certain covariates on the excess
mortality. Some covariates are cancer-specific whereas others are
variables that may influence the background mortality as well. The
latter should be taken into account in the background mortality to avoid
biases in estimating their effects on the excess mortality.
Unfortunately, the available life table might not include such variables
and, consequently, might provide inaccurate values of the background
mortality. We propose a model that uses multiplicative parameters to
correct potentially inaccurate background mortality. The model can be
seen as an extension of the frequently used Estève model because we
assume a Cox-type model for the excess mortality with a piecewise
constant baseline function and introduce additional parameters that
multiply the background mortality. The original and the extended model
are compared, first in a simulation study, then in an application to
colon cancer registry data.

A related software package can be found at a gitlab webpage or at
<https://CRAN.R-project.org/package=xhaz>.

## Installation

The most recent version of `xhaz` can be installed directly from the cran
repository using 

    install.packages("xhaz")

`xhaz` depends on the `stats`, `survival`, `optimParallel`, `optimx`,
`numDeriv`, `statmod`, `gtools` and `splines` packages which can be
installed directly from CRAN.

It also utilizes `survexp.fr`, the R package containing the french life
table. For example, to install `survexp.fr` follow the instructions
available at the [RStudio page on R and
survexp.fr](https://cran.r-project.org/package=survexp.fr).

First, install the R package via github

    devtools::install_github("rstudio/survexp.fr")

Then, when these other packages are installed, please load the `xhaz` R
package.

``` r
library(xhaz)
```

    ## Le chargement a nécessité le package : statmod

    ## Le chargement a nécessité le package : survival

### Fitting an classical excess hazard model with a piecewise constant baseline hazard

We illustrate the Esteve model using a simulated dataset from the
original Touraine et al. (2020) paper. This dataset is comprised of
2,000 patients with an information regarding their race as this
information can impact the patient background mortality. The US life
table can be used for the estimation of the model parameters.

``` r
data("simuData", package = "xhaz")

head(simuData)
```

    ##        age       agec  sex  race       date     time status time_year
    ## 1 50.52825 -1.3186092 male black 1990-10-21 72.00000      0  6.000000
    ## 2 62.50834 -0.4380647 male black 1990-07-17 72.00000      0  6.000000
    ## 3 64.49190 -0.2922719 male black 1990-10-18 20.15947      1  1.679956
    ## 4 48.23570 -1.4871131 male white 1990-12-06 13.78553      1  1.148794
    ## 5 31.71262 -2.7015697 male black 1990-01-20 72.00000      0  6.000000
    ## 6 31.11493 -2.7455005 male black 1990-01-17 18.73011      1  1.560842

``` r
dim(simuData)
```

    ## [1] 2000    8

``` r
levels(simuData$sex) <- c("male", "female")
interval <- c(0, 0.718, 1.351, 2.143, 3.601, 6)

fit.estv1 <- xhaz(formula = Surv(time_year, status) ~ agec + race,
                  data = simuData,
                  ratetable = survexp.us,
                  interval = interval,
                  rmap = list(age = 'age', sex = 'sex', year = 'date'),
                  baseline = c("constant"),
                  pophaz = "classic")
                  
fit.estv1
```

    ## Call:
    ## xhaz(formula = Surv(time_year, status) ~ agec + race, data = simuData, 
    ##     ratetable = survexp.us, rmap = list(age = "age", sex = "sex", 
    ##         year = "date"), baseline = c("constant"), pophaz = "classic", 
    ##     interval = interval)
    ## 
    ## 
    ##                coef se(coef) lower 0.95 upper 0.95     z Pr(>|z|)    
    ## agec         0.3340   0.0384     0.2588     0.4092  8.70  < 2e-16 ***
    ## racewhite   -0.5737   0.1422    -0.8524    -0.2949 -4.03  5.5e-05 ***
    ## [0-0.72[     0.1577   0.0124     0.1335     0.1820 12.74  < 2e-16 ***
    ## [0.72-1.35[  0.2317   0.0171     0.1981     0.2652 13.52  < 2e-16 ***
    ## [1.35-2.14[  0.2226   0.0168     0.1897     0.2556 13.24  < 2e-16 ***
    ## [2.14-3.6[   0.1489   0.0122     0.1250     0.1728 12.21  < 2e-16 ***
    ## [3.6-6[      0.1216   0.0107     0.1007     0.1426 11.37  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## 
    ## Excess hazard hazard ratio(s)
    ## (proportional effect variable(s) for exess hazard ratio(s))
    ##           exp(coef) lower 0.95 upper 0.95
    ## agec         1.3965     1.2953     1.5056
    ## racewhite    0.5635     0.4264     0.7446
    ## 
    ## number of observations: 2000;  number of events: 1375
    ## Likelihood ratio test: 103  on 7 degree(s) of freedom, p=0

### Fitting an excess hazard model with a piecewise constant baseline hazard with background mortality corrected with proportional effect for race variable

The new parameter to be added to `xhaz()` function is `add.rmap`: it
allows to specify the additional variable for the life table needed for
the estimation of the excess hazard parameters. This model concerns that
proposed by Touraine et al (2020).

``` r
fit.corrected1 <- xhaz(formula = Surv(time_year, status) ~ agec + race,
                       data = simuData,
                       ratetable = survexp.us,
                       interval = interval,
                       rmap = list(age = 'age', sex = 'sex', year = 'date'),
                       baseline = "constant", pophaz = "corrected",
                       add.rmap = "race")
                       
fit.corrected1
```

    ## Call:
    ## xhaz(formula = Surv(time_year, status) ~ agec + race, data = simuData, 
    ##     ratetable = survexp.us, rmap = list(age = "age", sex = "sex", 
    ##         year = "date"), baseline = "constant", pophaz = "corrected", 
    ##     add.rmap = "race", interval = interval)
    ## 
    ## 
    ##                     coef se(coef) lower 0.95 upper 0.95     z Pr(>|z|)    
    ## agec              0.2942   0.0798     0.1379     0.4506  3.69  0.00023 ***
    ## racewhite        -0.2082   0.1941    -0.5887     0.1723 -1.07  0.28000    
    ## [0-0.72[          0.1428   0.0205     0.1026     0.1829  6.97  3.2e-12 ***
    ## [0.72-1.35[       0.2143   0.0248     0.1657     0.2629  8.65  < 2e-16 ***
    ## [1.35-2.14[       0.2067   0.0264     0.1549     0.2585  7.83  5.0e-15 ***
    ## [2.14-3.6[        0.1324   0.0240     0.0853     0.1796  5.51  3.6e-08 ***
    ## [3.6-6[           0.1087   0.0238     0.0620     0.1554  4.56  5.1e-06 ***
    ## log(alpha.black)  0.3132   0.3068    -0.2880     0.9145  1.02  0.31000    
    ## log(alpha.white) -0.9172   1.0303    -2.9365     1.1021 -0.89  0.37000    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## 
    ## Excess hazard hazard ratio(s)
    ## (proportional effect variable(s) for exess hazard ratio(s))
    ##           exp(coef) lower 0.95 upper 0.95
    ## agec         1.3421     1.1478      1.569
    ## racewhite    0.8121     0.5551      1.188
    ## 
    ## and corrected scale parameters on population hazard 
    ##             exp(coef) lower 0.95 upper 0.95
    ## alpha.black    1.3679    0.74973      2.496
    ## alpha.white    0.3996    0.05305      3.010
    ## 
    ## number of observations: 2000;  number of events: 1375
    ## Likelihood ratio test: 14  on 9 degree(s) of freedom, p=0.124

### Fitting an excess hazard model with a piecewise constant baseline hazard with background mortality corrected with non proportional effect for race variable

The new parameter to be added to `xhaz` function is `add.rmap.cut`: it
furthermore allows to specify that the additional variable have a non
proportional effect on the background mortality. This excess hazard
model concerns that proposed by Mba et al (2020).

``` r
 # An additionnal cavariate (here race) missing in the life table is
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
```

    ## Call:
    ## xhaz(formula = Surv(time_year, status) ~ agec + race, data = simuData, 
    ##     ratetable = survexp.us, rmap = list(age = "age", sex = "sex", 
    ##         year = "date"), baseline = "constant", pophaz = "corrected", 
    ##     add.rmap = "race", add.rmap.cut = list(breakpoint = TRUE, 
    ##         cut = 75), interval = interval)
    ## 
    ## 
    ##                               coef se(coef) lower 0.95 upper 0.95     z
    ## agec                        0.2699   0.1020     0.0701     0.4697  2.65
    ## racewhite                  -0.3291   0.2983    -0.9138     0.2557 -1.10
    ## [0-0.72[                    0.1248   0.0249     0.0760     0.1735  5.02
    ## [0.72-1.35[                 0.1771   0.0286     0.1211     0.2331  6.20
    ## [1.35-2.14[                 0.1705   0.0302     0.1114     0.2297  5.65
    ## [2.14-3.6[                  0.1034   0.0293     0.0460     0.1607  3.53
    ## [3.6-6[                     0.0709   0.0267     0.0185     0.1233  2.65
    ## log(alpha.black_(30.9,75])  1.2886   0.1955     0.9054     1.6719  6.59
    ## log(alpha.black_(75,90.9])  0.6061   0.2359     0.1437     1.0685  2.57
    ## log(alpha.white_(30.9,75])  1.0259   0.4242     0.1945     1.8573  2.42
    ## log(alpha.white_(75,90.9]) -0.2266   0.4818    -1.1708     0.7176 -0.47
    ##                            Pr(>|z|)    
    ## agec                        0.00810 ** 
    ## racewhite                   0.27000    
    ## [0-0.72[                    5.3e-07 ***
    ## [0.72-1.35[                 5.7e-10 ***
    ## [1.35-2.14[                 1.6e-08 ***
    ## [2.14-3.6[                  0.00041 ***
    ## [3.6-6[                     0.00800 ** 
    ## log(alpha.black_(30.9,75])  4.4e-11 ***
    ## log(alpha.black_(75,90.9])  0.01000 ** 
    ## log(alpha.white_(30.9,75])  0.01600 *  
    ## log(alpha.white_(75,90.9])  0.64000    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## 
    ## Excess hazard hazard ratio(s)
    ## (proportional effect variable(s) for exess hazard ratio(s))
    ##           exp(coef) lower 0.95 upper 0.95
    ## agec         1.3098      1.073      1.600
    ## racewhite    0.7196      0.401      1.291
    ## 
    ## and corrected scale parameters on population hazard 
    ##  (non proportional correction using breakpoint approach)
    ## 
    ##                       exp(coef) lower 0.95 upper 0.95
    ## alpha.black_(30.9,75]    3.6278     2.4729      5.322
    ## alpha.black_(75,90.9]    1.8333     1.1546      2.911
    ## alpha.white_(30.9,75]    2.7896     1.2147      6.406
    ## alpha.white_(75,90.9]    0.7972     0.3101      2.049
    ## 
    ## 
    ## number of observations: 2239;  number of events: 1491
    ## Likelihood ratio test: 7.24  on 11 degree(s) of freedom, p=0.779

We can compare the output of these two models using AIC or BIC criteria.

``` r
AIC(fit.estv1)
```

    ## [1] 6767.122

``` r
AIC(fit.corrected1)
```

    ## [1] 6766.209

``` r
BIC(fit.estv1)
```

    ## [1] 6753.122

``` r
BIC(fit.corrected1)
```

    ## [1] 6748.209

A statistical comparison between two nested models can be performed with
a likelihood ratio test calculated by function `anova` method
implemented in `xhaz`.

As an example, say that we want to test whether we can drop all the
complex terms in the Mba model compared to the Touraine model.

As in `survival` package, We compare the two models using `anova()`,
i.e.,

``` r
anova(fit.corrected1, fit.corrected2)
```

    ## Assumption: Model 1 nested within Model 2
    ## 
    ## Likelihood ratio test
    ## Model 1: 
    ## Surv(time_year, status) ~ agec + race
    ## Model 2: 
    ## Surv(time_year, status) ~ agec + race
    ##      Model.df  loglik      df  Chisq Pr(>Chisq)    
    ## [1,]      9.0 -3374.1      NA     NA         NA    
    ## [2,]     11.0 -3528.7     2.0 154.56  < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Note that the user is responsible to supply appropriately nested excess
hazard models such that the LRT to be valid. The result suggests that we
could use the Mba model with non-proportional population hazards,
correcting for the life table with additional stratification on the
variable `race`.

### Plot of net survival and excess hazard for different models

One could be interested to the prediction of net survival and excess
hazard of the individual with the same characteristics as individual 1
on the the `simuData` (age 50.5, race black)

``` r
#only add Surv variables (time_year and status) to have them in the new.data. 
#They are not used for the prediction

simuData[1,]
```

    ##        age      agec  sex  race       date time status time_year
    ## 1 50.52825 -1.318609 male black 1990-10-21   72      0         6

``` r
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
```

![](README_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
par(old.par)
```

One could be interested to the prediction of marginal net survival and
marginal excess hazard of the individual with the same characteristics
as observed in the `simuData`.

``` r
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
```

![](README_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
par(old.par)
```

## License

GPL 3.0, for academic use.

## Acknowledgments

We are grateful to the members of the CENSUR Survival Group for their
helpful comments.
