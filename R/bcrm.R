## Bayesian CRM - extending original code of J. Jack Lee and Nan Chen,
## Department of Biostatistics,  the University of Texas M. D. Anderson Cancer Center
## Now using exact inference,  rjags,  R2WinBUGS or BRugs

### MJS 20/11/15

#' Bayesian Continual Reassessment Method for Phase I Dose-Escalation Trials
#'
#' Implements a wide variety of Bayesian CRM designs,  including 1-parameter,
#' 2-parameter and Escalation With Overdose Control (EWOC) designs. The program
#' can run interactively,  allowing the user to enter outcomes after each cohort
#' has been recruited,  or via simulation to assess operating characteristics.
#'
#' \code{bcrm} implements a Bayesian continual reassessment method (CRM)
#' (O'Quigley \emph{et al.},  1990); an adaptive design in which cohorts of
#' patients are sequentially recruited into a Phase I trial. A binary toxicity
#' outcome is assumed (e.g. Dose Limiting Toxicity / No Dose Limiting
#' Toxicity). The current cohort are given a dose "closest" to the specified
#' target toxicity level,  as estimated from the posterior distributions of
#' toxicity at each dose level from the patients thus far recruited. If
#' \code{pointest="mean"} then the posterior mean probability of toxicity is
#' used to choose the next dose. If \code{pointest="plugin"},  however,  the
#' posterior mean of the model parameter(s) is plugged-into the functional form
#' of the dose-toxicity model. To implement an EWOC design (Babb \emph{et al.},
#' 1998),  \code{pointest} should be a quantile,  \emph{q},  between 0 and 0.5.
#' The posterior distribution of the MTD (the dose in which the probability of
#' toxicity is equal to the target toxicity) is then calculated and the next
#' patient is given dose closest to the \emph{q}th quantile of the MTD
#' distribution.
#'
#' Alternatively,  escalation can be based on intervals of toxicity from the
#' posterior distribution using a loss function,  see Neuenschwander \emph{et
#' al.},  2008. To implement this approach,  the user should specify the
#' cutpoints of the toxicity intervals using \code{tox.cutpoints} and the
#' associated losses using \code{loss}.
#'
#' The possible choice of dose-toxicity model can be specified using \code{ff},
#' and includes the 1-parameter hyberbolic tangent,  logistic or power "working
#' models",  and the 2-parameter logistic model as follows: \describe{
#' \item{Hyperbolic
#' Tangent}{\deqn{p(Tox|d^*)=\left[(tanh(d^*)+1)/2\right]^\alpha}{p(Tox|d*)=[(tanh(d*)+1)/2]^\alpha}}
#' \item{Logistic (1-parameter)}{\deqn{p(Tox|d^*)=\frac{\exp(3+\alpha
#' d^*)}{1+\exp(3+\alpha d^*)}}{p(Tox|d*)=exp(3+\alpha d*)/(1+exp(3+\alpha
#' d*))}} \item{Power}{\deqn{p(Tox|d^*)={d^*}^\alpha}{p(Tox|d*)=d*^\alpha}}
#' \item{Logistic
#' (2-parameter)}{\deqn{p(Tox|d^*)=\frac{\exp(\log(\alpha_1)+\alpha_2
#' d^*)}{1+\exp(\log(\alpha_1)+\alpha_2 d^*)}}{
#' p(Tox|d*)=exp(log(\alpha_1)+\alpha_2 d*)/(1+exp(log(\alpha_1)+\alpha_2
#' d*))}} } where \eqn{\alpha>0} is the single positive-valued parameter for
#' the 1-parameter models,  and \eqn{\log(\alpha_1)}{log(\alpha_1)} and
#' \eqn{\alpha_2>0} are the intercept and slope parameters of the 2-parameter
#' model.
#'
#' The standardised doses,  \eqn{d^*}{d*},  are specified by the user using
#' \code{sdose},  or alternatively the prior probability of toxicity at each
#' dose level is specified using \code{p.tox0}. If the latter is used,  then the
#' standardised doses are calculated using the inverse of the functional form
#' and a plug-in estimate of the prior mean or median,  as specified in
#' \code{sdose.calculate},  as follows \deqn{d^* = f^{-1}(\code{p.tox0}, \alpha=
#' a)}{d* = f^{-1}(p.tox0, \alpha= a)} where \eqn{f^{-1}} is the the inverse of
#' the chosen functional form,  and the parameter(s) of the model are set equal
#' to \eqn{a},  either the prior mean or median of \eqn{\alpha}.
#'
#' Data that have already been accrued can be entered using the \code{data}
#' argument. A constrained CRM design can be implemented using
#' \code{constrain=TRUE},  in which case dose-skipping is prohibited (i.e. the
#' next cohort can only be dosed up to one dose level above the current
#' cohort). If a constrained model is used then the starting dose must be
#' specified using \code{start}. Alternatively,  if data have already been
#' accrued,  then the dose level of the last recruited patient determines the
#' constraint for the next patient.
#'
#' The prior is set using \code{prior.alpha}. For example
#' \code{prior.alpha=list(1, 1, 1)} specifies a Gamma prior with shape and scale
#' parameters both equal to one (\emph{i.e.} an Exponential(1) distribution),
#' whilst \code{prior.alpha=list(2, 0, 10)} specifies a Uniform(0, 10) prior.
#'
#' To specify a fixed maximum sample size of size \code{m} use
#' \code{stop=list(nmax=m)}. Alternatively,  the trial can stop after \code{m2}
#' patients have been treated at the current MTD estimate,  by setting
#' \code{stop=list(nmtd=m2)}. To stop the trial when the MTD estimate is within
#' a certain level of precision,  use \code{stop=list(precision=c(l, u))},  where
#' \code{l} and \code{u} are the lower and upper percentage points that the MTD
#' 95\% credible intervals for the risk of toxicity should lie within. Finally,
#' to prevent the trial stopping too early using these rules,  the argument
#' \code{stop=list(nmin=m3)} can be used to ensure the sample size is greater
#' than or equal to \code{m3}. Stopping rules can be used on their own or in
#' combination.
#'
#' The trial can be run interactively using \code{simulate=FALSE},  where the
#' user enters the outcomes for each new cohort,  or as a simulation study when
#' \code{simulate=TRUE}.
#'
#' The default calculations use exact methods (\code{method="exact"}) to
#' calculate the mean and quantiles for the posterior distributions. There are
#' three choices for MCMC calculations: \code{method="rjags"},
#' \code{method="BRugs"} or \code{method="R2WinBUGS"}. The first uses the JAGS
#' software,  the second uses OpenBUGS,  whilst the latter uses WinBUGS. To
#' implement these methods,  users require one or more of these packages to be
#' installed on their system.
#'
#' A simulated \code{bcrm} design can be compared with the standard 3+3
#' rule-based method,  see \code{\link{threep3}} for more details.
#'
#' @param stop A list of stopping rules for the trial. One or more of the
#' following options should be specified \describe{ \item{list("nmax")}{ The
#' maximum sample size of the trial} \item{list("nmtd")}{ The maximum number to
#' be treated at final maximum tolerated dose (MTD) estimate,  \emph{i.e.} if
#' the next recommended dose has already been administered to \code{nmtd}
#' patients then the trial will stop} \item{list("precision")}{ A vector of the
#' lower and upper percentage points that the MTD 95\% credible intervals for
#' the risk of toxicity should lie within} \item{list("nmin")}{ The minimum
#' sample size of the trial. To be used in conjunction with \code{nmtd} or
#' \code{precision}} }
#' @param data A named data frame giving information about dose and toxicity
#' from previously recruited patients. If missing,  then it is assumed that no
#' data have thus far been collected. Contains the following variables:
#' \describe{ \item{list("patient")}{ Recruited patient numbers,
#' \code{1, ..., n}} \item{list("dose")}{ Dose levels of recruited patients,
#' ranging from \code{1, ..., k}} \item{list("tox")}{ An indicator variable for
#' each patient (1=toxicity,  0=no toxicity)} }
#' @param p.tox0 A vector of length \code{k} listing the prior probabilities of
#' experiencing the outcome at each dose level \code{1, ...k}. The standardised
#' dose levels (CRM "skeleton") are formed from these probabilities using the
#' inverse of the functional form,  with a plug-in estimate for the prior mean
#' or median of alpha,  as specified in \code{ff},  \code{prior.alpha} and
#' \code{sdose.calculate}. Alternatively standardised dose levels can be given
#' directly using \code{sdose}.
#' @param sdose A vector of length \code{k} listing the standardised doses to
#' be used in the CRM model. Only required if \code{p.tox0} is missing.
#' @param dose Optional vector of length \code{k} of actual doses for plotting
#' purposes
#' @param ff A string indicating the functional form of the dose-response
#' curve. Options are \describe{ \item{ht}{ 1-parameter hyperbolic tangent}
#' \item{logit1}{ 1-parameter logistic} \item{power}{ 1-parameter power}
#' \item{logit2}{ 2-parameter logistic} }
#' @param prior.alpha A list of length 3 containing the distributional
#' information for the prior. The first element is a number from 1-4 specifying
#' the type of distribution. Options are \enumerate{ \item Gamma(a, b),  where
#' a=shape,  b=scale: mean=a*b,  variance=a*b*b \item Uniform(a, b),  where a=min,
#' b=max \item Lognormal(a, b),  where a=mean on the log scale,  b=variance on the
#' log scale \item Bivariate Lognormal(a, b),  where a=mean vector on the log
#' scale,  b=Variance-covariance matrix on the log scale. This prior should be
#' used only in conjunction with a two-parameter logistic model.  } The second
#' and third elements of the list are the parameters a and b,  respectively.
#' @param cohort The size of each cohort of patients that are sequentially
#' recruited to the trial. Defaults to 3
#' @param target.tox The target toxicity probability. Defaults to 1/3.
#' @param constrain Should a dose-skipping constraint be placed on the
#' escalation procedure,  as imposed by a modified CRM? Defaults to TRUE.
#' @param sdose.calculate What plug-in estimate of the prior alpha should be
#' used to calculate the standardised doses? Options are \code{"mean"}
#' (default) or \code{"median"}. Only required if \code{sdose} is missing.
#' @param pointest Which summary estimate of the posterior distribution should
#' be used to choose the next dose. Options are \code{"plugin"} (default) where
#' the posterior mean of the model parameter(s) is plugged into the function
#' form to obtain estimates of toxicity,  or \code{"mean"} where the posterior
#' mean probabilities of toxicity are directly used. Alternatively,  a number
#' between 0 and 1 can be specified representing the quantile of the maximum
#' tolerated dose (MTD) posterior distribution (e.g. 0.5 specifies the posteior
#' median). This produces an Escalation With Overdose Control (EWOC) design if
#' the quantile is less than 0.5 (see details). Currently,  EWOC designs must be
#' fit using MCMC methods.
#' @param tox.cutpoints A vector of cutpoints for toxicity intervals if these
#' are to be used to choose next dose. Defaults to NULL. For example
#' Underdosing [0, 0.2],  Target dosing (0.2,  0.35],  Excessive toxicity (0.35,
#' 0.60],  Unacceptable toxicity (0.60,  1.00] set
#' \code{tox.cutpoints=c(0.2, 0.35, 0.60)}.
#' @param loss A vector of length \code{length(tox.cutpoints)+1} specifying the
#' losses associated with each toxicity interval,  e.g. Underdosing = 1,  Target
#' dosing =0,  Excessive toxicity=1,  Unacceptable toxicity=2
#' @param start Dose level used at the beginning of the trial. Required if
#' \code{constrain=TRUE}.
#' @param simulate Should a simulation be conducted to assess operating
#' characteristics? Defaults to \code{TRUE}. If \code{FALSE},  a single CRM
#' trial is run interactively,  allowing the user to input outcomes after each
#' cohort is recruited.
#' @param nsims Number of simulations to perform if \code{simulate==TRUE}
#' (defaults to 1).
#' @param truep A vector of length k giving the true probabilities of the
#' outcome (toxicity) at each dose level \code{1, ..., k} in order to simulate
#' data. Only required if \code{simulate=TRUE}
#' @param threep3 Should operating characteristics of a standard 3+3 rule-based
#' design be calculated,  for comparison with \code{bcrm} design? Defaults to
#' \code{FALSE}. Only used in a simulation study,  i.e. when
#' \code{simulate=TRUE}.
#' @param method Optimisation method: options are \code{"exact"} (the default),
#' \code{"rjags"},  \code{"BRugs"},  or \code{"R2WinBUGS"}.
#' @param burnin.itr Number of burn-in iterations (default 2000).
#' @param production.itr Number of production iterations (default 2000).
#' @param bugs.directory Directory that contains the WinBUGS executable if
#' \code{method="R2WinBUGS"}. Defaults to "C:/Program Files/WinBUGS14/".
#' @param plot Should the dose-response curve be plotted after each cohort has
#' been entered? Defaults to FALSE.
#' @param seed Integer defining the state of the random number generator to
#' allow reproducibility of results. The default is to not specify a seed.
#' @param quietly Should the simulation number output be suppressed when
#' running bcrm?  Defaults to FALSE.
#' @param file File name where the dose-response plots are stored,  in a pdf
#' format. The program will ammend the current sample size to the end of the
#' file name.
#' @param N Final sample size (deprecated). To be replaced with \code{stop} in
#' future versions.
#' @param tox (Deprecated). A vector of length \code{k} listing the number of
#' patients who have experienced the outcome (toxicity) at each dose level
#' \code{1, ..., k}.
#' @param notox (Deprecated). A vector of length \code{k} listing the number of
#' patients who have not experienced the outcome (toxicity) at each dose level
#' \code{1, ..., k}.
#' @return \code{bcrm} returns an object of class "bcrm" or "bcrm.sim"; the
#' latter occuring when a simulation has been conducted (\code{simulate=TRUE}).
#' The function \code{\link{print}} (i.e. \code{\link{print.bcrm}} or
#' \code{\link{print.bcrm.sim}}) can be used to obtain summary information
#' about the design used,  the data observed,  current posterior estimates of
#' toxicity,  and the next recommended dose level.
#'
#' An object of class "bcrm" is a list with the following components:
#' \item{dose}{Range of doses} \item{sdose}{Standardised doses} \item{tox}{A
#' vector of length \code{k} listing the number of patients who have
#' experienced the outcome (toxicity) at each dose level \code{1, ..., k}}
#' \item{notox}{A vector of length \code{k} listing the number of patients who
#' have not experienced the outcome (toxicity) at each dose level
#' \code{1, ..., k}} \item{ndose}{A list of lists containing for each cohort the
#' components \code{ndose},  the dose level recommended for the next patient,
#' \code{est},  the estimated probabilitieis of toxicity using the chosen metric
#' (e.g. plugin,  mean,  quantile),  \code{mean},  the posterior mean probability
#' of toxicity at each dose,  \code{sd},  the posterior standard deviation for
#' probability of toxicity at each dose,  \code{quantiles},  the posterior
#' quantiles for probability of toxicity at each dose. This information is only
#' provided for cohorts recruited subsequent to any data given using \code{tox}
#' and \code{notox}. The first component relates to the prior information.}
#' \item{constrain}{Whether a constrained CRM design was used} \item{start}{The
#' starting dose for the latest run of the model if \code{constrain=TRUE}}
#' \item{target.tox}{The target toxicity level} \item{ff}{A number from 1-4
#' identifying the functional form,  1 = Hyperbolic tangent,  2 = 1-parameter
#' logistic,  3 = Power,  4 = 2-parameter logistic} \item{method}{The calculation
#' method used} \item{pointest}{The summary estimate used to choose the next
#' dose,  see \code{pointest}} \item{prior.alpha}{Information about the prior
#' used for \eqn{alpha},  see \code{prior.alpha}} \item{data}{A data frame with
#' variables `patient',  `dose' and `tox' listing the doses and outcomes of all
#' patients in the trial}
#'
#' An object of class "bcrm.sim" is a list of length \code{nsims}. Each
#' component is itself a list with components similar to those obtained from a
#' "bcrm" object. The print function,  \code{\link{print.bcrm.sim}} should be
#' used to obtain operating characteristics from the simulation.
#' @note Currently,  the re-parameterisation of the two-parameter model proposed
#' by (Babb \emph{et al.},  1998) is not implemented. Therefore,  users wishing
#' to implement an EWOC design should check whether their choice of prior for
#' the model parameter(s) translates to a sensible prior for the MTD
#' distribution before they implement the design. For example \preformatted{
#' prior.alpha <- list(1, 1, 1) ff <- "ht" target.tox <- 0.2
#' samples.alpha <- getprior(prior.alpha, 2000)
#' mtd <- find.x(ff, target.tox, alpha=samples.alpha) hist(mtd) }
#'
#' One-parameter models are designed as working models only,  and should not be
#' used with an escalation strategy based on intervals of the posterior
#' probabilities of toxicity.
#' @author Michael Sweeting \email{mjs212@@medschl.cam.ac.uk} (University of
#' Cambridge,  UK),  drawing on code originally developed by J. Jack Lee and Nan
#' Chen,  Department of Biostatistics,  the University of Texas M. D. Anderson
#' Cancer Center
#' @seealso \code{\link{print.bcrm}},  \code{\link{print.bcrm.sim}},
#' \code{\link{plot.bcrm}},  \code{\link{plot.bcrm.sim}},  \code{\link{threep3}}
#' @references Sweeting M.,  Mander A.,  Sabin T. \pkg{bcrm}: Bayesian Continual
#' Reassessment Method Designs for Phase I Dose-Finding Trials. \emph{Journal
#' of Statistical Software} (2013) 54: 1--26.
#' \url{http://www.jstatsoft.org/article/view/v054i13}
#'
#' O'Quigley J.,  Pepe M.,  Fisher L. Continual reassessment method: a practical
#' design for phase I clinical trials in cancer. \emph{Biometrics} (1990) 46:
#' 33--48.
#'
#' Babb J.,  Rogatko A.,  Zacks S. Cancer phase I clinical trials: efficient dose
#' escalation with overdose control. \emph{Statistics in Medicine} (1998) 17:
#' 1103--1120.
#'
#' Neuenschwander B.,  Branson M.,  Gsponer T. Critical aspects of the Bayesian
#' approach to phase I cancer trials. \emph{Statistics in Medicine} (2008) 27:
#' 2420--2439.
#' @examples
#'
#' ## Dose-escalation cancer trial example as described in Neuenschwander et al 2008.
#' ## Pre-defined doses
#' dose <- c(1, 2.5, 5, 10, 15, 20, 25, 30, 40, 50, 75, 100, 150, 200, 250)
#' ## Pre-specified probabilities of toxicity
#' ## [dose levels 11-15 not specified in the paper,  and are for illustration only]
#' p.tox0 <- c(0.010, 0.015, 0.020, 0.025, 0.030, 0.040, 0.050, 0.100, 0.170, 0.300, 0.400, 0.500, 0.650
#'   , 0.800, 0.900)
#' ## Data from the first 5 cohorts of 18 patients
#' data <- data.frame(patient=1:18, dose=rep(c(1:4, 7), c(3, 4, 5, 4, 2)), tox=rep(0:1, c(16, 2)))
#' ## Target toxicity level
#' target.tox <- 0.30
#'
#' ## A 1-parameter power model is used,  with standardised doses calculated using
#' ## the plug-in prior median
#' ## Prior for alpha is lognormal with mean 0 (on log scale)
#' ## and standard deviation 1.34 (on log scale)
#' ## The recommended dose for the next cohort if posterior mean is used
#' Power.LN.bcrm <- bcrm(stop=list(nmax=18), data=data, p.tox0=p.tox0, dose=dose
#'   , ff="power", prior.alpha=list(3, 0, 1.34^2), target.tox=target.tox, constrain=FALSE
#'   , sdose.calculate="median", pointest="mean")
#' print(Power.LN.bcrm)
#' plot(Power.LN.bcrm)
#'
#' ## Simulate 10 replicate trials of size 36 (cohort size 3) using this design
#' ## with constraint (i.e. no dose-skipping) and starting at lowest dose
#' ## True probabilities of toxicity are set to pre-specified probabilities (p.tox0)
#' Power.LN.bcrm.sim <- bcrm(stop=list(nmax=36), p.tox0=p.tox0, dose=dose, ff="power"
#'   , prior.alpha=list(3, 0, 1.34^2), target.tox=target.tox, constrain=TRUE
#'   , sdose.calculate="median", pointest="mean", start=1, simulate=TRUE, nsims=10, truep=p.tox0)
#' print(Power.LN.bcrm.sim)
#' plot(Power.LN.bcrm.sim)
#'
#' ## Comparing this CRM design with the standard 3+3 design
#' ## (only considering the first 12 dose levels)
#' \dontrun{
#' Power.LN.bcrm.compare.sim <- bcrm(stop=list(nmax=36), p.tox0=p.tox0[1:12], dose=dose[1:12]
#'   , ff="power", prior.alpha=list(3, 0, 1.34^2), target.tox=target.tox, constrain=TRUE
#'   , sdose.calculate="median", pointest="mean", start=1, simulate=TRUE, nsims=50
#'   , truep=p.tox0[1:12], threep3=TRUE)
#' print(Power.LN.bcrm.compare.sim, threep3=TRUE)
#' plot(Power.LN.bcrm.compare.sim, threep3=TRUE)
#' }
#'
#' ## A 2-parameter model,  using priors as specified in Neuenschwander et al 2008.
#' ## Posterior mean used to choose the next dose
#' ## Standardised doses using reference dose,  250mg
#' sdose <- log(dose/250)
#' ## Bivariate lognormal prior for two parameters
#' mu <- c(2.15, 0.52)
#' Sigma <- rbind(c(0.84^2, 0.134), c(0.134, 0.80^2))
#' ## Using rjags (requires JAGS to be installed)
#' TwoPLogistic.mean.bcrm <- bcrm(stop=list(nmax=18), data=data, sdose=sdose
#'   , dose=dose, ff="logit2", prior.alpha=list(4, mu, Sigma), target.tox=target.tox
#'   , constrain=FALSE, pointest="mean", method="rjags")
#' print(TwoPLogistic.mean.bcrm)
#' plot(TwoPLogistic.mean.bcrm)
#'
#' ## A 2-parameter model,  using an EWOC design with feasibility bound (MTD quantile)
#' ## of 0.25 to choose the next dose
#' ## Using rjags (requires JAGS to be installed)
#' \dontrun{
#' TwoPLogistic.EWOC0.25.bcrm <- bcrm(stop=list(nmax=18), data=data, sdose=sdose, dose=dose
#' 	  , ff="logit2", prior.alpha=list(4, mu, Sigma), target.tox=target.tox, constrain=FALSE
#' 	  , pointest=0.25, method="rjags")
#' print(TwoPLogistic.EWOC0.25.bcrm)
#' plot(TwoPLogistic.EWOC0.25.bcrm)
#' }
#'
#' ## A 2-parameter model,  using a loss function based on intervals of toxicity to choose
#' ## the next dose
#' ## Using rjags (requires JAGS to be installed)
#' \dontrun{
#' ## Toxicity cut-points
#' tox.cutpoints <- c(0.2, 0.35, 0.6)
#' ## Losses associated with toxicity intervals
#' ## [0, 0.2]=1,  (0.2, 0.35]=0,  (0.35, 0.6]=1,  (0.6, 1]=2
#' loss <- c(1, 0, 1, 2)
#' TwoPLogistic.tox.intervals.bcrm <- bcrm(stop=list(nmax=18), data=data, sdose=sdose
#'   , dose=dose, ff="logit2", prior.alpha=list(4, mu, Sigma), target.tox=target.tox
#'   , constrain=FALSE, tox.cutpoints=tox.cutpoints, loss=loss, method="rjags")
#' print(TwoPLogistic.tox.intervals.bcrm)
#' plot(TwoPLogistic.tox.intervals.bcrm)
#' ## Greater loss associated with overdosing and unacceptable toxicity
#' ## [0, 0.2]=1,  (0.2, 0.35]=0,  (0.35, 0.6]=2,  (0.6, 1]=4
#' loss2 <- c(1, 0, 2, 4)
#' TwoPLogistic.tox.intervals.2.bcrm <- bcrm(stop=list(nmax=18), data=data, sdose=sdose
#'   , dose=dose, ff="logit2", prior.alpha=list(4, mu, Sigma), target.tox=target.tox
#'   , constrain=FALSE, tox.cutpoints=tox.cutpoints, loss=loss2, method="rjags")
#' print(TwoPLogistic.tox.intervals.2.bcrm)
#' plot(TwoPLogistic.tox.intervals.2.bcrm)
#' }
#'
#'
#' @export bcrm
bcrm <- function(stop=list(nmax=NULL, nmtd=NULL, precision=NULL, nmin=NULL),
                 data=NULL, p.tox0=NULL, sdose=NULL, dose=NULL, ff=NULL, prior.alpha=NULL,
                 cohort=3, target.tox=NULL, constrain=TRUE, sdose.calculate="mean",
                 pointest="plugin", tox.cutpoints=NULL, loss=NULL, start=NULL,
                 simulate=FALSE, nsims=1, truep=NULL, threep3=FALSE,
                 method="exact", burnin.itr=2000, production.itr=2000,
                 bugs.directory="c:/Program Files/WinBUGS14/", plot=FALSE,
                 seed=NULL,  quietly=FALSE,  file=NULL, N=NULL, tox=NULL, notox=NULL, ...){

  wh <- check.bcrm.args(N, stop, tox, notox, p.tox0, sdose, sdose.calculate,
                        pointest, method, tox.cutpoints, simulate, truep,
                        ff, constrain, start, data, loss, prior.alpha)

  if(!is.null(tox.cutpoints)){
    pointest <- NULL
  }
  if(simulate){
    plot <- FALSE
    results <- subset.results <- list()
  }

  # Set seed if specified
  if(!is.null(seed)){
    set.seed(seed)
  }

  if(is.null(sdose)){
    sdose <- get.bcrm.sdose(prior.alpha, sdose.calculate, ff, p.tox0)
  }

  # Checking that length of truep in simulation study is same length as sdose
  if(simulate & length(truep) != length(sdose)) stop("Length of truep must be the same as length of sdose or p.tox0")
  # Checking that length of dose (if specified) is same length as sdose
  if(length(dose)>0 & length(dose)!=length(sdose)) stop("Length of dose must be the same as length of sdose or p.tox0")

  # Check data contains the correct variables
  if(!is.null(data)){
    wh <- check.bcrm.data(data, sdose, start)
    data <- data[order(data$patient), ]
    start <- as.numeric(data$dose[1])
  }

  ## Allowing access to fast exact computation if method="exact" & simulate=TRUE &
  ## stopping rules do not depend on posterior quantiles
  if(method == "exact" & simulate & is.null(stop$precision)){
    method <- "exact.sim"
  }
  ## If method=="exact" and a two-parameter model is fitted,  only relevant
  ## escalation posterior quantities are calculated (apart from trial end)
  if(method == "exact" & ff == "logit2"){
    method <- "exact.sim"
    if (plot){
      plot <- FALSE
      warning("Plot function not available for exact computation of 2-parameter model")
    }
  }
  k <- length(sdose)
  sim <- 1

  while(sim <= nsims){
    if(simulate){
      sub.start <- 100*floor((sim-1)/100)
      if(sub.start > length(results)){
        results <- c(results, subset.results)
        subset.results <- list()
      }
    }
    ## Prior
    if(is.null(data)){
      new.tox <-  rep(0, k)
      new.notox <- rep(0, k)
      current <- start-1
      alpha <- if(sim == 1) switch(method
      , rjags = getprior(prior.alpha,  10000)
      , BRugs = getprior(prior.alpha,  10000)
      , R2WinBUGS = getprior(prior.alpha, 10000)
      , exact = Posterior.exact(new.tox, new.notox, sdose, ff, prior.alpha)
      , exact.sim = Posterior.exact.sim(new.tox, new.notox, sdose, ff, prior.alpha, pointest)
      )
    } else {
      new.tox <- as.numeric(xtabs(tox~factor(dose, levels=1:k), data=data))
      new.notox <- as.numeric(xtabs((1-tox)~factor(dose, levels=1:k), data=data))
      current <- as.numeric(data$dose[dim(data)[1]])
      alpha <- if (sim == 1) switch(method
      , rjags = Posterior.rjags(new.tox, new.notox, sdose, ff, prior.alpha, burnin.itr, production.itr)
      , BRugs = Posterior.BRugs(new.tox, new.notox, sdose, ff, prior.alpha, burnin.itr, production.itr)
      , R2WinBUGS = Posterior.R2WinBUGS(new.tox, new.notox, sdose, ff, prior.alpha, burnin.itr, production.itr, bugs.directory)
      , exact = Posterior.exact(new.tox, new.notox, sdose, ff, prior.alpha)
      , exact.sim = Posterior.exact.sim(new.tox, new.notox, sdose, ff, prior.alpha, pointest)
      )
    }
    ncurrent <- sum(new.tox+new.notox)
    if(is.null(data)){
      newdata <- data.frame(patient=NULL, dose=NULL, tox=NULL)
    } else{
      newdata <- data
    }
    found <- FALSE
    match <- 1
    if (sim == 1){
      prior.ndose <- switch(method
        , rjags = nextdose(alpha, sdose, ff, target.tox, constrain, pointest, current, tox.cutpoints, loss)
        , BRugs = nextdose(alpha, sdose, ff, target.tox, constrain, pointest, current, tox.cutpoints, loss)
        , R2WinBUGS = nextdose(alpha, sdose, ff, target.tox, constrain, pointest, current, tox.cutpoints, loss)
        , exact = nextdose.exact(alpha, sdose, ff, target.tox, constrain, pointest, current)
        , exact.sim = nextdose.exact.sim(alpha, sdose, ff, target.tox, constrain, pointest, current)
        )
    }
    ndose <- prior.ndose

    if(!simulate){
      results <- list(dose=dose, sdose=sdose, tox=new.tox, notox=new.notox, ndose=list(ndose),
                      constrain=constrain, start=start, target.tox=target.tox, ff=ff,
                      method=method, pointest=pointest, tox.cutpoints=tox.cutpoints, loss=loss,
                      prior.alpha=prior.alpha, data=data)
      class(results) <- "bcrm"
    } else {
      subset.results[[sim-sub.start]] <- list(dose=dose, sdose=sdose, tox=new.tox, notox=new.notox,
                                              ndose=list(ndose), constrain=constrain, start=start,
                                              target.tox=target.tox, ff=ff, method=method,
                                              pointest=pointest, tox.cutpoints=tox.cutpoints, loss=loss,
                                              prior.alpha=prior.alpha, truep=truep)
    }

    if (plot){
      plot(results, file)
    }

    stopped <- stop.check(stop, ncurrent, ndose, new.tox, new.notox, simulate)

    while (!stopped){
      current <- ndose[[1]]
      ncurrent <- ncurrent+cohort
      if(!simulate){
        interact <- crm.interactive(new.tox, new.notox, ncurrent, cohort, ndose, dose)
        if(interact$bk==TRUE){
          results$data <- newdata
          return(results)
        }
        y <- interact$y
        new.tox <- interact$tox
        new.notox <- interact$notox
      } else {
        y <- rbinom(cohort, 1, truep[ndose[[1]]])
        new.notox[ndose[[1]]] <- new.notox[ndose[[1]]]+(cohort-sum(y))
        new.tox[ndose[[1]]] <- new.tox[ndose[[1]]]+sum(y)
      }
      currentdata <- data.frame(patient=(ncurrent-cohort+1):ncurrent, dose=rep(current, cohort), tox=y)
      newdata <- rbind(newdata, currentdata)

      if(simulate & match<length(results) & method!="exact.sim"){
        repeat{
          match.tox <- all(xtabs(tox~factor(dose, levels=1:k), data=results[[match]]$data[1:dim(newdata)[1], ])==xtabs(tox~factor(dose, levels=1:k), data=newdata))
          match.notox <- all(xtabs((1-tox)~factor(dose, levels=1:k), data=results[[match]]$data[1:dim(newdata)[1], ])==xtabs((1-tox)~factor(dose, levels=1:k), data=newdata))
          match.current <- results[[match]]$data$dose[dim(newdata)[1]]==newdata$dose[dim(newdata)[1]]
          m <- match.tox & match.notox & match.current
    # alternative  m <- all(results[[match]]$data[1:dim(newdata)[1], ] == newdata)
          if(m){
            ndose <- results[[match]]$ndose[[length(subset.results[[sim-sub.start]]$ndose)+1]]
            found <- TRUE
            break
          } else if(match!=(length(results)-1)){
            match <- match+1
          } else {
            match <- match+1
            found <- FALSE
            break
          }
        }
      }
      if(!found){
        alpha <- switch(method
          , rjags=Posterior.rjags(new.tox,  new.notox,  sdose,  ff,  prior.alpha,  burnin.itr,  production.itr)
          , BRugs=Posterior.BRugs(new.tox,  new.notox,  sdose,  ff,  prior.alpha,  burnin.itr,  production.itr)
          , R2WinBUGS=Posterior.R2WinBUGS(new.tox,  new.notox,  sdose,  ff,  prior.alpha,  burnin.itr,  production.itr, bugs.directory)
          , exact=Posterior.exact(new.tox, new.notox, sdose, ff, prior.alpha)
          , exact.sim=Posterior.exact.sim(new.tox, new.notox, sdose, ff, prior.alpha, pointest)
          )
        ndose <- switch(method
          , rjags=nextdose(alpha, sdose, ff, target.tox, constrain, pointest, current, tox.cutpoints, loss)
          , BRugs=nextdose(alpha, sdose, ff, target.tox, constrain, pointest, current, tox.cutpoints, loss)
          , R2WinBUGS=nextdose(alpha, sdose, ff, target.tox, constrain, pointest, current, tox.cutpoints, loss)
          , exact=nextdose.exact(alpha, sdose, ff, target.tox, constrain, pointest, current)
          , exact.sim=nextdose.exact.sim(alpha, sdose, ff, target.tox, constrain, pointest, current)
          )
      }
      if(!simulate){
        results$tox <- new.tox
        results$notox <- new.notox
        results$ndose[[length(results$ndose)+1]] <- ndose
        results$data <- newdata
      } else{
        subset.results[[sim-sub.start]]$tox <- new.tox
        subset.results[[sim-sub.start]]$notox <- new.notox
        subset.results[[sim-sub.start]]$ndose[[length(subset.results[[sim-sub.start]]$ndose)+1]] <- ndose
        subset.results[[sim-sub.start]]$data <- newdata      }
      if(plot){
        plot(results, file)
      }
      stopped <- stop.check(stop, ncurrent, ndose, new.tox, new.notox, simulate)
    } # Close while (!stopped)

    if(simulate & !quietly){
      message(sim)
    }
  sim <- sim+1
  } # Close while (sim < nsims)

  if(simulate){
    results <- c(results, subset.results)
    class(results) <- "bcrm.sim"
  }
  if(simulate & threep3){
    if(length(truep)>9){
      cat("\n Warning: Calculation of all 3+3 designs may take a long time,  continue?  ")
      yn  <-  readline()
                if (yn!="y" & yn=="Y") return(results)
          }
    cat("\n Calculating operating characteristics of a standard 3+3 trial for comparison... \n")
    results[[1]]$threep3 <- threep3(truep, start=start, dose=dose)
  }
  return(results)
}
