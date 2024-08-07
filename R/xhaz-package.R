#' @title Excess Hazard Modelling Considering Inappropriate Mortality Rates
#'
#' @aliases xhaz-package
#'
#' @description Contains functions to fit excess hazard models, with or without
#' proportional population hazards assumption. The baseline excess hazard could be a
#' piecewise constant function or a B-splines. When B-splines is choosen for
#' the baseline excess hazard, the user can specify some covariates which have
#' a time-dependent effect (using "bsplines") on the baseline excess hazard.
#' The user can also specify if the framework corresponds to the classical
#' excess hazard modeling, i.e. assuming that the expected mortality of studied
#' individuals is appropriate. He can also consider two other frameworks: first,
#' the expected mortality available in the life table is not accurate and
#' requires taking into account an additional variable in the life table by allowing the
#' latter acts on the general population morality with a proportional effect.
#' This approach is presented by Touraine et al. (2020) <doi:10.1177/0962280218823234>.
#' The user can also fit a model that relax the proportional expected hazards assumption
#' considered in the latter excess hazard model. This extension was proposed by
#' Mba et al. (2020) <doi:10.1186/s12874-020-01139-z> allows non-proportional
#' effect of the additional variable on the general population mortality;
#' second, there is a non-comparability source of bias in terms of expected mortality
#' of selected individuals in a non-population-based studies such as clinical trials.
#' The related excess hazard model correcting this source of bias is presented in
#' Goungounga et al. (2019) <doi:10.1186/s12874-019-0747-3>. The optimization process
#' in these presented models uses the maximum likelihood method through the routine
#' \code{optim} or an internal function of the \code{xhaz-package}.
#' Finally, it is now possible to fit a rescaled excess hazard regression model using
#' various flexible parametrisations proposed by the mexhaz R package.
#' Moreover, the package enables the implementation of a rescaled random effects
#' excess hazard regression model, which allows for the examination of heterogeneity
#' between clusters (i.e., random effects at the cluster level) impacts on excess hazard,
#' as well as the investigation of potentially inappropriate mortality rates. In other words,
#' the latter concerns situations where the assumption of comparability between other-cause
#' mortality rates in the cohort under study and observed mortality rates in the general
#' population may be invalid, which could introduce bias into the estimation of excess
#' hazard model parameters (Goungounga et al. (2023) <doi:10.1002/bimj.202100210>).
#'
#' @details
#' \tabular{ll}{
#' Package: \tab xhaz\cr
#' Type: \tab Package\cr
#' Version: \tab 2.0.2\cr
#' Date: \tab 2024-04-25\cr
#' License: \tab GPL-3\cr
#' }
#'
#' @references Goungounga JA, Touraine C, Graff\'eo N, Giorgi R;
#' CENSUR working survival group. Correcting for misclassification
#' and selection effects in estimating net survival in clinical trials.
#' BMC Med Res Methodol. 2019 May 16;19(1):104.
#' doi: 10.1186/s12874-019-0747-3. PMID: 31096911; PMCID: PMC6524224.
#' (\href{https://pubmed.ncbi.nlm.nih.gov/31096911/}{PubMed})
#'
#' Touraine C, Graff\'eo N, Giorgi R; CENSUR working survival group.
#' More accurate cancer-related excess mortality through correcting
#' background mortality for extra variables.
#' Stat Methods Med Res. 2020 Jan;29(1):122-136.
#' doi: 10.1177/0962280218823234. Epub 2019 Jan 23. PMID: 30674229.
#' (\href{https://pubmed.ncbi.nlm.nih.gov/30674229/}{PubMed})
#'
#' Mba RD, Goungounga JA, Graff\'eo N, Giorgi R; CENSUR working survival group.
#' Correcting inaccurate background mortality in excess hazard models
#' through breakpoints. BMC Med Res Methodol. 2020 Oct 29;20(1):268.
#' doi: 10.1186/s12874-020-01139-z. PMID: 33121436; PMCID: PMC7596976.
#' (\href{https://pubmed.ncbi.nlm.nih.gov/33121436/}{PubMed})
#'
#'
#' Goungounga, JA, Graff\'eo N, Charvat H, Giorgi R. “Correcting for
#' heterogeneity and non-comparability bias in multicenter clinical trials
#' with a rescaled random-effect excess hazard model.” Biometrical journal.
#' Biometrische Zeitschrift vol. 65,4 (2023): e2100210.
#' doi:10.1002/bimj.202100210.PMID: 36890623;
#' (\href{https://pubmed.ncbi.nlm.nih.gov/36890623/}{PubMed})
#'
#'
#' @examples
#'
#' \donttest{
#' library("numDeriv")
#' library("survexp.fr")
#' library("splines")
#' data("simuData", "dataCancer", package = "xhaz")
#' # load the data sets 'simuData' and 'dataCancer'.
#'
#'
#'# Esteve et al. model: baseline excess hazard is a piecewise function
#'#                      linear and proportional effects for the covariates on
#'#                      baseline excess hazard.
#'
#'
#' fit.estv1 <- xhaz(formula = Surv(time_year, status) ~ agec + race,
#'                   data = simuData,
#'                   ratetable = survexp.us,
#'                   interval = c(0, NA, NA, NA, NA, NA, max(simuData$time_year)),
#'                   rmap = list(age = 'age', sex = 'sex', year = 'date'),
#'                   baseline = "constant",
#'                   pophaz = "classic")
#'
#'
#' fit.estv1
#'
#'
#' # Touraine et al. model: baseline excess hazard is a piecewise function
#' #                        with a linear and proportional effects for the
#' #                        covariates on the baseline excess hazard.
#' # An additionnal cavariate (here race) missing in the life table is
#' # considered by the model.
#'
#'
#' fit.corrected1 <- xhaz(formula = Surv(time_year, status) ~ agec + race,
#'                        data = simuData,
#'                        ratetable = survexp.us,
#'                        interval = c(0, NA, NA, NA, NA, NA,
#'                                     max(simuData$time_year)),
#'                        rmap = list(age = 'age', sex = 'sex', year = 'date'),
#'                        baseline = "constant", pophaz = "corrected",
#'                        add.rmap = "race")
#'
#'
#'
#' fit.corrected1
#'
#'
#'
#'  # An additionnal cavariate (here race) missing in the life table is
#'  # considered by the model with a breakpoint at 75 years
#'
#'  fit.corrected2 <- xhaz(formula = Surv(time_year, status) ~ agec + race,
#'                         data = simuData, ratetable = survexp.us,
#'                         interval = c(0, NA, NA, NA, NA, 6),
#'                         rmap = list(age = 'age', sex = 'sex', year = 'date'),
#'                         baseline = "constant", pophaz = "corrected",
#'                         add.rmap = "race",
#'                         add.rmap.cut = list(breakpoint = TRUE, cut = 75))
#'
#'
#'  fit.corrected2
#'
#' data("breast")
#' #load the data sets 'breast'.
#'
#'  # Flexible mexhaz model: baseline excess hazard with cubic B-splines
#'  # assumption on the life table available :
#'  # other cause mortality in the cohort is comparable to the mortality
#'  # observed in the general population with the same characteristics.
#'
#' # The life table to be used is survexp.us. Note that SEX is coded 2 instead of female in survexp.us.
#' breast$sexe <- "female"
#'
#' fit.haz <- exphaz(
#'                   formula = Surv(temps, statut) ~ 1,
#'                   data = breast, ratetable = survexp.us,
#'                   only_ehazard = FALSE,
#'                   rmap = list(age = 'age', sex = 'sexe', year = 'date'))
#'
#' breast$expected <- fit.haz$ehazard
#' breast$expectedCum <- fit.haz$ehazardInt
#'
#' mod.bs <- mexhazLT(formula = Surv(temps, statut) ~ agecr + armt,
#'                   data = breast,
#'                   ratetable = survexp.us, degree = 3,
#'                   knots=quantile(breast[breast$statut==1,]$temps, probs=c(1:2/3)),
#'                   expected = "expected",expectedCum = "expectedCum",
#'                   base = "exp.bs", pophaz = "classic")
#'
#' mod.bs
#'
#'
#'  # Flexible mexhaz model: baseline excess hazard with cubic B-splines
#'  # assumption on the life table available :
#'  # other cause mortality in the cohort is different to the mortality
#'  # observed in the general population with the same characteristics.
#'
#' mod.bs2 <- mexhazLT(formula = Surv(temps, statut) ~ agecr + armt,
#'                   data = breast, degree = 3,
#'                   knots=quantile(breast[breast$statut==1,]$temps, probs=c(1:2/3)),
#'                   expected = "expected",expectedCum = "expectedCum",
#'                   base = "exp.bs", pophaz = "rescaled")
#'
#' mod.bs2
#'
#'
#'  # Flexible mexhaz model with a random effects at cluster level:
#'  # baseline excess hazard with cubic B-splines
#'  # assumption on the life table used :
#'  # other cause mortality in the cohort is different to the mortality
#'  # observed in the general population with the same characteristics.
#'
#' mod.bs3 <- mexhazLT(formula = Surv(temps, statut) ~ agecr + armt,
#'                   data = breast, degree = 3,
#'                   knots=quantile(breast[breast$statut==1,]$temps, probs=c(1:2/3)),
#'                   expected = "expected",expectedCum = "expectedCum",
#'                   base = "exp.bs", pophaz = "rescaled", random = "hosp")
#'
#' mod.bs3
#'
#'
#' }
#'
#'
#' @keywords internal
"_PACKAGE"


#> [1] "_PACKAGE"
