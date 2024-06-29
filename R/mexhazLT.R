#' @title mexhazLT function
#'
#' @description
#' Extends excess hazard models from the mexhaz R`-`package to allow rescaling (Goungounga et al. (2019) <doi:10.1186/s12874-019-0747-3>) of the background mortality in the presence or absence of multilevel data (Goungounga et al. (2023) <doi: 10.1002/bimj.202100210>).
#' It allows for different shapes of the baseline hazard, the ability to include time`-`dependent effects of variable(s), and a random effect at the cluster level.
#'
#' @param formula a formula object of the  function with the response on the left
#'  of a \code{~} operator and the terms on the right. The response must be a
#'  survival object as returned by the \code{Surv} function (time in first and status
#'  in second).
#'
#' @note `time` is OBLIGATORY in YEARS.
#'
#' @param data a data frame in which the variables named in the formula are to be interpreted.
#'
#' @param expected name of the variable (must be given in quotes) representing
#' the population instantaneous hazard.
#'
#' @param expectedCum name of the variable (must be given in quotes) representing
#' the population cumulative hazard.
#'
#' @param base functional form that should be used to model the baseline hazard.
#' Selection can be made between the following options: \code{"weibull"} for a Weibull
#' hazard, \code{"exp.bs"} for a hazard described by the exponential of a B`-`spline
#' (only B`-`splines of degree 1, 2 or 3 are accepted), \code{"exp.ns"} for a hazard
#' described by the exponential of a restricted cubic spline (also called 'natural spline'), \code{"pw.cst"} for a piecewise constant hazard. By default, base="weibull" as
#' in mexhaz R`-`package.
#'
#' @param pophaz specifies two possible arguments in character: classic and
#'  rescaled. If \code{pophaz = "classic"} is chosen, it fits the models that do
#'  not require the background mortality  to be rescaled and assumes that the
#'  comparability assumption holds; if \code{pophaz = "rescaled"} is chosen, it
#'  fits the models that require that require the background mortality to be rescaled.
#'
#' @param degree if \code{base="exp.bs"}, degree represents the degree of the B`-`spline used.
#' Only integer values between 1 and 3 are accepted, and 3 is the default.
#'
#' @param knots if \code{base="exp.bs"} or \code{"exp.ns"}, knots is the vector
#' of interior knots of the spline. If \code{base="pw.cst"}, knots is the vector
#' defining the endpoints of the time intervals on which the hazard is assumed
#' to be constant. By default, \code{knots=NULL} (that is, it produces a B`-`spline
#' with no interior knots if base="exp.bs", a linear B`-`spline with no interior
#' knots if base="exp.ns", or a constant hazard over the whole follow`-`up period
#' if \code{base="pw.cst"}).
#'
#' @param bound a vector of two numerical values corresponding to the boundary
#' knots of the spline functions. If \code{base="exp.bs"} or \code{base="exp.ns"},
#' computation of the B-spline basis requires that boundary knots be given.
#' The bound argument allows the user to specify these boundary knots.
#' If \code{base="exp.bs"}, the interval defined by the boundary knots must at
#' least include the interval \code{c(0,max(time))} (otherwise, there could be
#'  problems with ill`-`conditioned bases). If \code{base="exp.ns"},
#'
#' @param n.gleg corresponds to the number of quadrature nodes to be specified as in \code{mexhaz}.
#'
#' @param init vector of initial values as in \code{mexhaz}.
#'
#' @param random name of the variable to be entered as a random effect (must be
#' given between quotes), representing the cluster membership. As in \code{mexhaz}
#' \code{random=NULL} means that the function fits a fixed effects model.
#'
#' @param n.aghq corresponds to the number of quadrature points to be specified
#' as in \code{mexhaz} for the estimation of the cluster`-`specific marginal
#' likelihoods by adaptative Gauss`-`Hermite quadrature.
#'
#' @param fnoptim name of the R optimisation procedure used to maximise the
#' likelihood. Selection can be made between "nlm" (by default) and "optim".
#' Note: if \code{exactGradHess=TRUE}, this argument will be ignored
#' (fnoptim will be set automatically to \code{"nlm"}).
#'
#'
#' @param verbose integer parameter representing the frequency at which the
#' current state of the optimisation process is displayed. If verbose=0 (default),
#' nothing is displayed.
#'
#' @param method if fnoptim="optim", method represents the optimisation method to
#' be used by optim. By default, \code{method="Nelder-Mead"}. This parameter is not
#' used if \code{fnoptim="nlm"}.
#'
#' @param iterlim if \code{fnoptim="nlm"}, iterlim represents the maximum number of
#' iterations before the nlm optimisation procedure is terminated. By default,
#' iterlim is set to 10000. This parameter is not used if \code{fnoptim="optim"}
#' (in this case, the maximum number of iterations must be given as part of a
#' list of control parameters via the control argument: see the help page of optim
#' for further details).
#'
#' @param numHess logical value allowing the user to choose between the Hessian
#' returned by the optimization algorithm (default) or the Hessian estimated by
#' the hessian function from the \code{numDeriv} package.
#'
#' @param print.level	 his argument is only used if \code{fnoptim="nlm"}. It determines
#' the level of printing during the optimisation process. The default value
#' (for the mexhaz function) is set to \code{'1'} which means that details on the
#' initial and final step of the optimisation procedure are printed (see the
#' help page of nlm for further details).
#'
#' @param exactGradHess logical value allowing the user to decide whether
#' maximisation of the likelihood should be based on the analytic gradient and
#' Hessian computed internally (default, corresponding to \code{exactGradHess=TRUE}).
#'
#' @param gradtol this argument is only used if \code{fnoptim="nlm"}.
#' It corresponds to the tolerance at which the scaled gradient is considered
#' close enough to zero to terminate the algorithm. The default value depends on
#' the value of the argument \code{exactGradHess}.
#'
#' @param testInit this argument is used only when \code{exactGradHess=TRUE} and
#' when the model is not an excess hazard random effect model. It instructs
#' the mexhaz function to try several vectors of initial values in case
#' optimization was not successful with the default (or user-defined) initial
#' values. Because optimization based on the analytical gradient and Hessian
#' is usually fast, this simple and empirical procedure proves useful to
#' increase the probability of convergence in cases when it is difficult
#' to specify appropriate initial values.
#'
#' @param keep.data logical argument determining whether the dataset should be
#' kept in the object returned by the function: this can be useful in certain
#' contexts (e.g., to calculate cluster`-`specific posterior predictions from a
#' random intercept model) but might create unnecessarily voluminous objects.
#' The default value is set to \code{FALSE}.
#'
#' @param ... other parameters used with the \code{mexhazLT} function
#'
#' @keywords mexhazLT
#'
#' @aliases mexhazAlpha
#'
#' @return An object of class \code{mexhaz}, \code{xhaz} or \code{mexhazLT}.
#' This object is a list containing the following components:
#'
#'
#' \item{dataset}{name of the dataset used to fit the model.}
#'
#' \item{call}{function call on which the model is based.}
#'
#' \item{formula}{formula part of the call.}
#'
#' \item{withAlpha}{logical value indicating whether the model corresponds to a class of models correcting for life tables.}
#'
#' \item{expected}{name of the variable corresponding to the population hazard.}
#'
#' \item{expectedCum}{name of the variable corresponding to the cumulative population hazard.}
#'
#' \item{xlevels}{information concerning the levels of the categorical variables used in the model.}
#'
#' \item{n.obs.tot}{total number of observations in the dataset.}
#'
#' \item{n.obs}{number of observations used to fit the model (after exclusion of missing values).}
#'
#' \item{n.events}{number of events (after exclusion of missing values).}
#'
#' \item{n.clust}{number of clusters.}
#'
#' \item{n.time.0}{number of observations for which the observed follow-up time was equal to 0 (only for right censored type data).}
#'
#' \item{base}{function used to model the baseline hazard.}
#'
#' \item{max.time}{maximal observed time in the dataset.}
#'
#' \item{boundary.knots}{vector of boundary values used to define the B`-`spline (or natural spline) bases.}
#'
#' \item{degree}{degree of the B`-`spline used to model the logarithm of the baseline hazard.}
#'
#' \item{knots}{vector of interior knots used to define the B`-`spline (or natural spline) bases.}
#'
#' \item{names.ph}{names of the covariables with a proportional effect.}
#'
#' \item{random}{name of the variable defining cluster membership (set to NA in the case of a purely fixed effects model).}
#'
#' \item{init}{a vector containing the initial values of the parameters.}
#'
#' \item{coefficients}{a vector containing the parameter estimates.}
#'
#' \item{std.errors}{a vector containing the standard errors of the parameter estimates.}
#'
#' \item{vcov}{the variance-covariance matrix of the estimated parameters.}
#'
#' \item{gradient}{the gradient of the log`-`likelihood function evaluated at the estimated parameters.}
#'
#' \item{hessian}{the Hessian of the log`-`likelihood function evaluated at the estimated parameters.}
#'
#' \item{mu.hat}{a data.frame containing the estimated cluster`-`specific random effects (shrinkage estimators).}
#'
#' \item{var.mu.hat}{the covariance matrix of the cluster`-`specific shrinkage estimators.}
#'
#' \item{vcov.fix.mu.hat}{a matrix containing the covariances between the fixed effect and the cluster`-`specific shrinkage estimators. More specifically, the i`-`th line of the matrix represents the covariances between the shrinkage estimator of the i`-`th cluster and the fixed effect estimates. This matrix is used by the function \code{predict.mexhaz} to make cluster`-`specific predictions.}
#'
#' \item{data}{original dataset used to fit the model (if \code{keep.data} was set to \code{TRUE}).}
#'
#' \item{n.par}{number of estimated parameters.}
#'
#' \item{n.gleg}{number of Gauss`-`Legendre quadrature points used to calculate the cumulative (excess) hazard (only relevant if a B-spline of degree 2 or 3 or a cubic restricted spline was used to model the logarithm of the baseline hazard).}
#'
#' \item{n.aghq}{number of adaptive Gauss`-`Hermite quadrature points used to calculate the cluster-specific marginal likelihoods (only relevant if a multi-level model is fitted).}
#'
#' \item{fnoptim}{name of the R optimisation procedure used to maximise the likelihood.}
#'
#' \item{method}{optimisation method used by optim.}
#'
#' \item{code}{code (integer) indicating the status of the optimisation process (this code has a different meaning for nlm and for optim: see their respective help page for details).}
#'
#' \item{loglik}{value of the log`-`likelihood at the end of the optimisation procedure. Note that this is different to that calculated in mexhaz as the cumulative expected hazard cannot be removed from the log`-`likelihood.}
#'
#' \item{iter}{number of iterations used in the optimisation process.}
#'
#' \item{eval}{number of evaluations used in the optimisation process.}
#'
#' \item{time.elapsed}{total time required to reach convergence.}
#'
#' @author
#' Juste Goungounga, Hadrien Charvat, Nathalie Graffeo, Roch Giorgi
#'
#' @references
#' Goungounga JA, Touraine C, Graff\'eo N, Giorgi R; CENSUR working
#'  survival group. Correcting for misclassification and selection effects in
#'  estimating net survival in clinical trials. BMC Med Res Methodol. 2019 May
#'  16;19(1):104. doi: 10.1186/s12874-019-0747-3. PMID: 31096911; PMCID:
#'  PMC6524224. (\href{https://pubmed.ncbi.nlm.nih.gov/31096911/}{PubMed})
#'
#'  Goungounga, JA, Graff\'eo N, Charvat H, Giorgi R. “Correcting for
#'  heterogeneity and non-comparability bias in multicenter clinical trials
#'  with a rescaled random-effect excess hazard model.” Biometrical journal.
#'  Biometrische Zeitschrift vol. 65,4 (2023): e2100210.
#'  doi:10.1002/bimj.202100210.PMID: 36890623;
#'  (\href{https://pubmed.ncbi.nlm.nih.gov/36890623/}{PubMed})
#'
#'
#'
#'
#' @examples
#' \donttest{
#' library("numDeriv")
#' library("survexp.fr")
#' library("splines")
#' library("statmod")
#' data("breast")
#' # load the data sets 'breast'.
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
#'}
#'
#'
#'
#' @import survival
#' @import stats
#' @import statmod
#' @import splines
#' @import survexp.fr
#' @import mexhaz
#'
#' @export
mexhazLT <- function (formula, data, expected = "expected", expectedCum = "expectedCum",
                      pophaz = "classic",
                      base = c("weibull", "exp.bs", "exp.ns", "pw.cst"),
          degree = 3, knots = NULL, bound = NULL,
          n.gleg = 20, init = NULL, random = NULL, n.aghq = 10, fnoptim = c("nlm", "optim"), verbose = 0,
          method = "Nelder-Mead", iterlim = 10000, numHess = FALSE,
          print.level = 1, exactGradHess = TRUE,
          gradtol = ifelse(exactGradHess,1e-08, 1e-06), testInit = TRUE,
          keep.data = FALSE, ...)
{

    resLT <- formulaLT(formula = formula, data = data, expected = expected,
                     expectedCum = expectedCum, pophaz = pophaz,
                     base = base, degree = degree, knots = knots, bound = bound,
                     n.gleg = n.gleg, init = init, random = random, n.aghq = n.aghq,
                     verbose = verbose, iterlim = iterlim,
                     print.level = print.level, gradtol = gradtol, testInit = testInit,
                     keep.data = keep.data, ...)
    resLT$call <- match.call()
    resLT$recurrent <- FALSE
    class(resLT)[2:3] <- c("xhaz", "mexhazLT")
    return(resLT)
}


