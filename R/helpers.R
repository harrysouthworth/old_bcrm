check.bcrm.data <- function(){
  if(any(!(c("patient", "dose", "tox") %in% names(data)))){
    stop("data must have variables named 'patient',  'dose' and 'tox'")
  }
  if(any(data$patient != 1:dim(data)[1])){
    stop("'patient' variable in data must be an ascending vector of positive integers")
  }
  if(any(!(data$tox %in% c(0, 1)))){
    stop("'tox' variable in data must be a vector of zeros (no toxicity) and ones (toxicity)")
  }
  if(any(!(data$dose %in% 1:length(sdose)))){
    stop(paste("'dose' variable in data must contain the dose levels (1 to ", length(sdose), ")", sep=""))
  }
  if(!is.null(start)){
    warning("start no longer needs to be specified if data is given; using last recruited patient as current dose")
  }
  invisible()
}

get.bcrm.sdose <- function(){
  alpha.prior.plug <- if (prior.alpha[[1]] == 1){
    ifelse(sdose.calculate == "mean", prior.alpha[[2]] * prior.alpha[[3]], median(getprior(prior.alpha,  10000)))
  } else if(prior.alpha[[1]] == 2){
    0.5*(prior.alpha[[2]]+prior.alpha[[3]])
  } else if(prior.alpha[[1]] == 3){
    ifelse(sdose.calculate=="mean", exp(prior.alpha[[2]] + prior.alpha[[3]] / 2), exp(prior.alpha[[2]]))
  } else if(prior.alpha[[1]] == 4){
    if(sdose.calculate == "mean"){
      exp(prior.alpha[[2]] + diag(prior.alpha[[3]])/2)
    } else {
        exp(prior.alpha[[2]])
    }
  }

  find.x(ff, p.tox0, alpha=alpha.prior.plug)
}

check.bcrm.args <- function(){
  # Checks of argument inputs
  if(missing(N) & is.null(stop$nmax) & is.null(stop$nmtd) & is.null(stop$precision)){
    stop("At least one stopping rule must be provided using the stop argument")
  }
  if(!missing(N)){
    stop$nmax <- N
    warning("N is deprecated and users should now use the stop argument to specify the maximum sample size")
  }
  if(!missing(tox) | !missing(notox)){
    stop("tox and nontox arguments are deprecated and users should now use the data argument to specify previous data,  see ?bcrm")
  }
  if(!(length(stop$precision) %in% c(0, 2))){
    stop("stop$precision must be a vector of length two")
  }
  if(!is.null(stop$nmax) & !is.null(stop$nmin)) {
    if(stop$nmin>stop$nmax) stop("stop$nmin must be less than stop$nmax")
  }
  if(missing(p.tox0) & missing(sdose)){
    stop("Either p.tox0 or sdose must be specified")
  }
  if(!missing(p.tox0) & !missing(sdose)){
    stop("Only one of p.tox0 and sdose must be specified")
  }
  if(sdose.calculate!="mean" & sdose.calculate!="median"){
    stop("sdose.calculate must be either `mean' or `median'")
  }
  if((is.character(pointest) & pointest!="mean" & pointest!="plugin") | is.numeric(pointest) & (pointest<0 | pointest>1)){
    stop("pointest must be either `plugin',  `mean' or an EWOC feasibility quantile between 0 and 1")
  }
  if(is.numeric(pointest) & method=="exact"){
    stop("EWOC design must be fitted using MCMC methods")
  }
  if(!is.null(tox.cutpoints) & method=="exact"){
    stop("Escalation based on toxicity intervals must be fit using MCMC. Please specify either method=`rjags',  method='BRugs' or method='R2WinBUGS'")
  }

  if(simulate & is.null(truep)){
    stop("truep must be specified if simulating data")
  }

  if(!(method %in% c("exact", "rjags", "BRugs", "R2WinBUGS"))){
    stop("method must be either `exact',  `rjags',  `BRugs' or `R2WinBUGS'")
  }

  ## Check to see if ff is one of "ht", "logit1", "power", "logit2"
  if((!ff %in% c("ht", "logit1", "power", "logit2"))){
    stop("ff must be one of `ht',  `logit1',  `power' or `logit2'")
  }

  if(ff=="logit2" & method=="exact"){
    warning("Exact method slow for 2-parameter model,  suggest using rjags (MCMC)")
  }
  if(constrain & is.null(start) & is.null(data)){
    stop("A starting dose level must be specified using `start' if constrain==TRUE")
  }
  if((!is.null(tox.cutpoints) & is.null(loss)) | (is.null(tox.cutpoints) & !is.null(loss))){
    stop("Both tox.cutpoints and loss must be specified to conduct escalation based on toxicity intervals")
  }
  if(!is.null(tox.cutpoints) & length(loss)!=length(tox.cutpoints)+1){
    stop("The number of losses must be one more than the number of cutpoints")
  }

  if(!is.null(tox.cutpoints) & ff!="logit2"){
    warning("One-parameter models are designed as working models only,  and should not be used with an escalation strategy based on intervals of the posterior probabilities of toxicity")
  }

  if(ff=="logit2" & (length(prior.alpha[[2]])<2 | length(prior.alpha[[3]])<2)){
    stop("second and third components of `prior.alpha' must be vectors of size 2")
  }

  invisible()
}

# ----------------------------------------------------------------------
#     User inputs toxicity info from the treatment at a dose level:
#   Must be either 0 and 1
# ----------------------------------------------------------------------
get.answer  <-  function() {
  repeat {
    ans  <-  as.numeric(readline())
    if (is.na(ans)) next
    if (ans != floor(ans)) next
    if (ans>=0 & ans<=1) return(ans)
  }
}

# ----------------------------------------------------------------------
#     ask the user to input an integer number <= n.
#     the number is checked to belong to [-1, n] and also to be
#     an integer.  `ENTER' returns to the caller with no action,
#
#     n - biggest number to be accepted
#    ncurrent - patient no. for current patient
#     cohort - cohort size
# ----------------------------------------------------------------------
get.dose.level  <-  function( n , ncurrent, cohort) {
  repeat {
    cat("\n\n ENTER DOSE LEVEL BETWEEN 1 AND ",  n, " FOR PATIENTS ", ncurrent-cohort+1, " to ", ncurrent)
    cat("\n (`RETURN' TO ACCEPT RECOMMENDATION,  0 TO EXIT AND RETURN CURRENT RESULTS)  ")

    ans  <-  readline()
    if ( nchar(ans)==0 ) return( -2 )
    ans  <-  as.integer(ans)
    if (is.na(ans)) next
    if ( -1<=ans && ans<=n ) return( ans )
  }
}

# ----------------------------------------------------------------------
#     ask the user to input a dose from those given
#     `ENTER' returns to the caller with no action,
#
#     dose  - dose labels
#    ncurrent - patient no. for current patient
#     cohort - cohort size
# ----------------------------------------------------------------------
get.dose  <-  function( dose , ncurrent, cohort) {
  repeat {
    cat("\n\n ENTER DOSE FOR PATIENTS ", ncurrent-cohort+1, " to ", ncurrent)
    cat("\n POSSIBLE CHOICES ARE ", dose)
    cat("\n (`RETURN' TO ACCEPT RECOMMENDATION,  0 TO EXIT AND RETURN CURRENT RESULTS)  ")

    ans  <-  readline()
    if ( nchar(ans)==0 ) return( -2 )
    if (is.na(ans)) next
    if ( ans %in% dose | ans==0) return( ans )
  }
}

# ----------------------------------------------------------------------
# generates a vector of values and prior distribution
#
#     prior.alpha --> list containing the information for prior
#               [[1]] - the prior distribution type:
#                     1 - gamma: mean=a*b; var=a*b*b
#                     2 - uniform: a+b*unif(0, 1)
#           3 - lognormal: lnorm(a, b),  where mean=a,  var=b
#           4 - log Multivariate normal(a, b),  where a=mean vector,  b=Variance-covariance matrix
#               [[2]] - a: first parameter of the prior distribution
#               [[3]] - b: second parameter of the prior distribution
#
# ----------------------------------------------------------------------


#' Samples from the specified prior distribution.
#'
#' A sample of specified size is obtained from the prior distribution.
#'
#' A vector of size \code{n} is returned from the specified prior distribution.
#'
#' @param prior.alpha A list of length 3 containing the distributional
#' information for the prior. The first element is a number from 1-4 specifying
#' the type of distribution. Options are \enumerate{ \item Gamma(a, b),  where
#' a=shape,  b=scale: mean=a*b,  variance=a*b*b \item Uniform(a, b),  where a=min,
#' b=max \item Lognormal(a, b),  where a=mean on the log scale,  b=variance on the
#' log scale \item Bivariate Lognormal(a, b),  where a=mean vector on the log
#' scale,  b=Variance-covariance matrix on the log scale. This prior should be
#' used only in conjunction with a two-parameter logistic model.  } The second
#' and third elements of the list are the parameters a and b,  respectively.
#' @param n The number of samples.
#' @author Michael Sweeting \email{mjs212@@medschl.cam.ac.uk} (University of
#' Cambridge,  UK),  drawing on code originally developed by J. Jack Lee and Nan
#' Chen,  Department of Biostatistics,  the University of Texas M. D. Anderson
#' Cancer Center
#' @seealso \code{\link{bcrm}},  \code{\link{find.x}}
#' @references Sweeting M.,  Mander A.,  Sabin T. \pkg{bcrm}: Bayesian Continual
#' Reassessment Method Designs for Phase I Dose-Finding Trials. \emph{Journal
#' of Statistical Software} (2013) 54: 1--26.
#' \url{http://www.jstatsoft.org/article/view/v054i13}
#' @examples
#'
#' prior.alpha <- list(1, 1, 1)
#' samples.alpha <- getprior(prior.alpha, 2000)
#' hist(samples.alpha)
#'
#' @export getprior
getprior  <-  function(prior.alpha,  n) {
  type  <-  prior.alpha[[1]]
  a     <-  prior.alpha[[2]]
  b     <-  prior.alpha[[3]]
  if ( type == 1 ) {
    prior <- rgamma(n, a) * b
  }
  else if ( type == 2 ) {
    prior <- sapply(1:length(a), function(i){ runif(n, a[i], b[i])} )
  }
  else if (type == 3) {
    prior <- rlnorm(n, a, sqrt(b))
  }
  else if (type == 4) {
    log.prior <- rmvnorm(n, a, b)
    prior <- exp(log.prior)
  }
  return (prior)
}
