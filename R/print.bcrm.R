#' Print information regarding a trial conducted using the Bayesian continuous
#' reassessment method
#'
#' Print method for a trial or series of trials conducted using a
#' \code{\link{bcrm}} model.
#'
#' If a single trial is conducted,  then the \code{\link{print}} function
#' currently produces summary information about the design used,  the data
#' observed,  current posterior estimates of toxicity,  and the next recommended
#' dose level.  If a simulation study is conducted,  then the following
#' operating characteristics are printed: \describe{ \item{Experimentation
#' proportion}{Proportion of patients recruited to each dose,  and to each true
#' region of toxicity,  across the simulated trials} \item{Recommendation
#' proportion}{Proportion of trials that recommend each of the dose levels as
#' the final maximum tolerated dose (i.e. with toxicity "closest" to the target
#' toxicity level),  and the associated regions of true toxicity for the
#' recommended MTDs} } If \code{trajectories = TRUE} then the dose level
#' administered and outcome observed are returned as matrices for every patient
#' (column) in every simulation (row).  If \code{threep3 = TRUE} then the
#' operating characteristics of the standard 3+3 design are displayed alongside
#' those of the \code{bcrm} design (see \code{\link{threep3}} for more
#' details).
#'
#' @aliases print.bcrm print.bcrm.sim
#' @param x An object of class "bcrm" or "bcrm.sim" as returned by
#' \code{\link{bcrm}}
#' @param tox.cutpoints An optional argument passed to \code{print.bcrm.sim}
#' specifying the cutpoints of toxicity for which the operating characteristics
#' are to be categorised. Defaults to \code{seq(from=0, to=1, by=0.2)}
#' @param trajectories Should the individual simulation dose and outcome
#' trajectories be returned? Defaults to FALSE.
#' @param threep3 Should operating characteristics of a standard 3+3 rule-based
#' design be displayed alongside those from the \code{bcrm} design? Defaults to
#' \code{FALSE}.
#' @param ... Further arguments passed to or from other methods
#' @return The following two components are returned from
#' \code{print.bcrm.sim}: \item{exp}{A matrix with number of rows equal to the
#' number of doses,  and number of columns equal to the number of simulations.
#' Gives the experimentation proportions for each dose within each simulation.}
#' \item{rec}{A vector with length equal to the number of simulations,  giving
#' the recommended MTD for each simulation.}
#' @author Michael Sweeting \email{mjs212@@medschl.cam.ac.uk} (University of
#' Cambridge,  UK)
#' @seealso \code{\link{bcrm}},  \code{\link{threep3}}
#' @references Sweeting M.,  Mander A.,  Sabin T. \pkg{bcrm}: Bayesian Continual
#' Reassessment Method Designs for Phase I Dose-Finding Trials. \emph{Journal
#' of Statistical Software} (2013) 54: 1--26.
#' \url{http://www.jstatsoft.org/article/view/v054i13}
print.bcrm <- function(x, ...){
  cat(" Estimation method: ", x$method, "\n")
  ff.txt <- switch(x$ff
                   , ht="Hyperbolic Tangent"
                   , logit1="1-parameter logistic"
                   , power="1-parameter power"
                   , logit2="Two-parameter logistic")

  cat("\n Model: ", ff.txt, "\n")

  pa.txt <- switch(x$prior.alpha[[1]]
                   , "1"=paste0("Gamma(Shape:", x$prior.alpha[[2]], ", Scale:", x$prior.alpha[[3]], ")", sep="")
                   , "2"=paste0("Uniform(", x$prior.alpha[[2]], ", ", x$prior.alpha[[3]], ")", sep="")
                   , "3"=paste0("Lognormal(Mean:", x$prior.alpha[[2]], ", Variance:", x$prior.alpha[[3]], ")", sep="")
                   , "4"=paste0("Log Multivariate Normal"))
  cat("\n Prior: ", pa.txt, "\n")
  if(x$prior.alpha[[1]]==4){
    cat("Mean Vector: \n")
    print(x$prior.alpha[[2]])
    cat("\nVariance-Covariance Matrix: \n")
    print(x$prior.alpha[[3]])
  }
  tab1 <- x$sdose
  names(tab1) <- x$dose
  cat("\n Standardised doses (skeleton): \n")
  print(tab1)

  if(x$constrain) {
    if(is.null(x$dose)){
      cat("\n Modified (constrained) CRM used,  starting dose level: ", x$start, "\n")
    } else {
      cat("\n Modified (constrained) CRM used,  starting dose: ", x$dose[x$start], "\n")
    }
  } else { cat("\n Unmodified (unconstrained) CRM used \n") }

  if(!is.null(x$loss)){
    cat("\n Loss function given intervals of toxicity used to select next dose.")
    tab.lf <- x$loss
    names(tab.lf) <- levels(cut(0, breaks=c(0, x$tox.cutpoints, 1)))
    cat("\n Loss function: \n")
    print(tab.lf)
  } else if(x$pointest=="plugin"){
    cat("\n Plug-in estimate of probability of toxicity used to select next dose \n")
  } else if(x$pointest=="mean"){
    cat("\n Posterior mean estimate of probability of toxicity used to select next dose \n")
  } else {
    cat("\n", 100*x$pointest, "percentile of (standardised) MTD distribution used to select next dose")
    cat("\n", 100*x$pointest, "percentile is:", x$ndose[[length(x$ndose)]]$target, "\n")
  }

  tab <- rbind(x$tox+x$notox, x$tox)
  rownames(tab) <- c("n", "Toxicities")
  colnames(tab) <- x$dose
  names(dimnames(tab)) <- c("", "Doses")
  cat("\n Toxicities observed: \n")
  print(tab)

  if(x$method %in% c("exact", "exact.sim") & x$ff=="logit2"){
    tab2a <- signif(rbind(x$ndose[[length(x$ndose)]]$est), 3)
    rownames(tab2a) <- c("Point estimate")
  } else {
    tab2a <- signif(rbind(x$ndose[[length(x$ndose)]]$mean, x$ndose[[length(x$ndose)]]$sd, x$ndose[[length(x$ndose)]]$quantiles["50%", ]), 3)
    rownames(tab2a) <- c("Mean", "SD", "Median")
    tab2b <- signif(x$ndose[[length(x$ndose)]]$quantiles, 3)
    colnames(tab2b) <- x$dose
    names(dimnames(tab2b)) <- c("Quantiles", "Doses")
  }
  colnames(tab2a) <- x$dose
  names(dimnames(tab2a)) <- c("", "Doses")
  cat("\n Posterior estimates of toxicity: \n")
  print(tab2a)
  if(!(x$method %in% c("exact", "exact.sim") & x$ff=="logit2")){
    print(tab2b)
  }
  if(!is.null(x$loss)){
    tab3 <- rbind(x$ndose[[length(x$ndose)]]$est)
    colnames(tab3) <- x$dose
    cat("\n Posterior expected loss at each dose: \n")
    print(tab3)
    tab4 <- x$ndose[[length(x$ndose)]]$probs
    colnames(tab4) <- x$dose
    rownames(tab4) <- levels(cut(0, breaks=c(0, x$tox.cutpoints, 1)))
    cat("\n Posterior probability of dose being in each toxicity interval")
    names(dimnames(tab4)) <- c("Toxicity intervals", "Doses")
    print(tab4)
  } else if(x$pointest=="plugin"){
    tab3 <- signif(rbind(x$ndose[[length(x$ndose)]]$est), 3)
    colnames(tab3) <- x$dose
    cat("\n Plug-in estimates of toxicity: \n")
    print(tab3)
  }
  if(is.null(x$dose)){
    cat("\n Next recommended dose level: ", x$ndose[[length(x$ndose)]]$ndose, "\n")
  } else {
    cat("\n Next recommended dose: ", x$dose[x$ndose[[length(x$ndose)]]$ndose], "\n")
  }
}

#-----------------------------------------------------------------------
#    Print function for an object of class bcrm.sim
# -----------------------------------
print.bcrm.sim <- function(x, tox.cutpoints=NULL, trajectories=FALSE, threep3=FALSE, ...){
  if(trajectories){
    ## sample size
    n <- sapply(x, function(i){dim(i$data)[1]})
    dose.mat <- outcome.mat <- matrix(NA, nrow=length(x), ncol=max(n))
    for(i in 1:length(x)){
      dose.mat[i, 1:n[i]] <- x[[i]]$data$dose
      outcome.mat[i, 1:n[i]] <- x[[i]]$data$tox
    }
    return(list(doses=dose.mat, outcomes=outcome.mat))
  } else {
    # average sample size
    n.average <- mean(sapply(x, function(i){dim(i$data)[1]}))
    n.min <- min(sapply(x, function(i){dim(i$data)[1]}))
    n.max <- max(sapply(x, function(i){dim(i$data)[1]}))
    if(n.min==n.max){
      tab0 <- cbind(n.average)
      rownames(tab0) <- "Sample size"
      colnames(tab0) <- ""
    } else {
      tab0 <- cbind(n.average, n.min, n.max)
      rownames(tab0) <- "Sample size"
      colnames(tab0) <- c("Mean", "Minimum", "Maximum")
    }
    exp <- sapply(x, function(i){(i$tox+i$notox)/sum(i$tox+i$notox)})
    exp.tab <- apply(exp, 1, mean)
    rec <- sapply(x, function(i){i$ndose[[length(i$ndose)]]$ndose})
    rec.tab <- prop.table(table(factor(rec, levels=1:length(x[[1]]$tox))))
    tab <- signif(rbind(exp.tab, rec.tab), 3)
    rownames(tab) <- c("Experimentation proportion", "Recommendation proportion")
    dose <- if(is.null(x[[1]]$dose)){1:length(x[[1]]$truep)} else {x[[1]]$dose}
    colnames(tab) <- dose
    names(dimnames(tab)) <- c("", "Doses")
    if(is.null(tox.cutpoints)){
      tox.cutpoints <- seq(0, 1, by=0.2)
    } else {
      tox.cutpoints <- unique(c(0, tox.cutpoints, 1))
    }
    exp.tox <- prop.table(table(cut(unlist(sapply(x, function(i){rep(i$truep, (i$tox+i$notox))}, simplify=FALSE)), tox.cutpoints, include.lowest=T)))
    rec.tox <- prop.table(table(cut(sapply(x, function(i){i$truep[i$ndose[[length(i$ndose)]]$ndose]}), tox.cutpoints, include.lowest=T)))
    tab2 <- signif(rbind(exp.tox, rec.tox), 3)
    rownames(tab2) <- c("Experimentation proportion", "Recommendation proportion")
    names(dimnames(tab2)) <- c("", "Probability of DLT")
    cat("Operating characteristics based on ", length(x), " simulations: \n \n")
    print(tab0)
    cat("\n")
    print(tab)
    cat("\n")
    print(tab2)
    if(threep3 & is.null(x[[1]]$threep3)){
      cat("\n Calculating 3+3 operating characteristics....\n")
      x[[1]]$threep3 <- threep3(x[[1]]$truep, x[[1]]$start)
    }
    if(threep3){
      cat("\n\n******************** 3+3 operating characteristics *****************\n")
      print.threep3(x[[1]]$threep3, tox.cutpoints=tox.cutpoints, dose=dose)
    }
  }
  invisible(list(rec=rec, exp=exp))
}
