#' Plot the operating characteristics from a standard 3+3 trial
#'
#' Plots of the operating characteristics obtained from a standard 3+3 trial,
#' using \code{\link{threep3}}
#'
#' This function plots the sample size distribution,  the experimentation
#' distribution,  the recommended dose distribution and the percentage of
#' subjects who experience the toxicity outcome (dose-limiting toxicity) for
#' the standard 3+3 trial.
#'
#' @param x An object of class "threep3",  as returned by \code{\link{threep3}}.
#' @param file File name where the operating characteristic plot is stored,  in
#' a pdf format.
#' @param ... Further arguments passed to or from other methods
#' @author Michael Sweeting \email{mjs212@@medschl.cam.ac.uk} (University of
#' Cambridge,  UK)
#' @seealso \code{\link{threep3}}
#' @references Sweeting M.,  Mander A.,  Sabin T. \pkg{bcrm}: Bayesian Continual
#' Reassessment Method Designs for Phase I Dose-Finding Trials. \emph{Journal
#' of Statistical Software} (2013) 54: 1--26.
#' \url{http://www.jstatsoft.org/article/view/v054i13}
plot.threep3 <- function(x, file=NULL, ...){
  dose <- if(is.null(x$dose)) 1:length(x$truep) else x$dose
  dose.label <- if(is.null(x$dose)) "Dose levels" else "Dose"

  # sample size
  n.threep3 <- x$ssize
  df.n.threep3 <- data.frame(n=n.threep3, weight=x$prob)
  a <- ggplot()+stat_bin(aes(x=n, y=100*..density.., weight=weight), data=df.n.threep3, binwidth=1)+xlab("Sample size")+ylab("Percent")+ggtitle("Sample size")

  # experimentation
  exp.threep3 <- rep(dose, 10000*x$exp)
  df.exp.threep3 <- data.frame(exp=as.factor(exp.threep3))
  b <- ggplot()+geom_bar(aes(x=exp, y=..count../100), data=df.exp.threep3)+xlab(dose.label)+ylab("Percent")+ggtitle("Experimentation")

  # recommendation
  rec.threep3 <- dose[x$mtd]
  df.rec.threep3 <- data.frame(rec=factor(rec.threep3), weight=x$prob[x$mtd!=0])
  c <- ggplot()+geom_bar(aes(x=rec, y=100*..count.., weight=weight), data=df.rec.threep3)+xlab(dose.label)+ylab("Percent")+ggtitle("Recommendation")

  # observed DLTs
  obs.threep3 <- 100*x$dlt.no/x$ssize
  bw <- max(diff(range(obs.threep3[x$prob>0]))/30, 1)
  df.obs.threep3 <- data.frame(obs=obs.threep3, weight=x$prob, bw=bw)
  df.obs.threep3 <- subset(df.obs.threep3, df.obs.threep3$weight>0)
  d <- ggplot()+geom_histogram(aes(bw=bw, x=obs, y=100*..density..*bw, weight=weight), data=df.obs.threep3, binwidth=bw)+xlab("Percentage of subjects with DLTs")+ylab("Percent")+ggtitle("DLTs")

  if(!is.null(file))
    pdf(paste(file, ".pdf", sep=""), ...)
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(2, 2)))
  vplayout <- function(x, y)  viewport(layout.pos.row=x, layout.pos.col=y)
  if(!is.null(a)) print(a, vp=vplayout(1, 1))
  print(b, vp=vplayout(1, 2))
  print(c, vp=vplayout(2, 1))
  print(d, vp=vplayout(2, 2))
  if(!is.null(file))
    dev.off()
}


#-----------------------------------------------------------------------
#    Print function for an object of class bcrm
# -----------------------------------


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
                   , "1"=paste("Gamma( Shape:", x$prior.alpha[[2]], ",  Scale:", x$prior.alpha[[3]], ")", sep="")
                   , "2"=paste("Uniform(", x$prior.alpha[[2]], ",  ", x$prior.alpha[[3]], ")", sep="")
                   , "3"=paste("Lognormal( Mean:", x$prior.alpha[[2]], ",  Variance:", x$prior.alpha[[3]], ")", sep="")
                   , "4"=paste("Log Multivariate Normal"))
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

#-----------------------------------------------------------------------
#    Print function for an object of class threep3
# -----------------------------------


#' Print information regarding the operating characteristics of a standard 3+3
#' design
#'
#' Print method for a 3+3 design specified using a \code{\link{threep3}}.
#'
#' The following operating characteristics are printed for the standard 3+3
#' design: \describe{ \item{Sample size}{Mean,  minimum and maximum sample size
#' of the design} \item{Experimentation proportion}{Proportion of patients
#' recruited to each dose,  and to each true region of toxicity,  on average}
#' \item{Recommendation proportion}{Proportion of 3+3 trials that would
#' recommend each of the dose levels as the final maximum tolerated dose (see
#' \code{\link{threep3}} for definition of the MTD),  and the associated regions
#' of true toxicity for the recommended MTDs} \item{Average number of
#' patients}{The average number of patients dosed at each level} \item{Average
#' number of DLTs}{The average number of DLTs seen at each level} }
#'
#' @param x An object of class "threep3" as returned by \code{\link{threep3}}
#' @param tox.cutpoints An optional argument passed to \code{print.threep3}
#' specifying the cutpoints of toxicity for which the operating characteristics
#' are to be categorised. Defaults to \code{seq(from=0, to=1, by=0.2)}
#' @param dose Optional vector of length \code{k} of actual doses for
#' presentation purposes
#' @param ... Further arguments passed to or from other methods
#' @author Michael Sweeting \email{mjs212@@medschl.cam.ac.uk} (University of
#' Cambridge,  UK)
#' @seealso \code{\link{threep3}}
#' @references Sweeting M.,  Mander A.,  Sabin T. \pkg{bcrm}: Bayesian Continual
#' Reassessment Method Designs for Phase I Dose-Finding Trials. \emph{Journal
#' of Statistical Software} (2013) 54: 1--26.
#' \url{http://www.jstatsoft.org/article/view/v054i13}
print.threep3 <- function(x, tox.cutpoints=NULL, dose=NULL, ...){
  if(is.null(dose)){
    dose <- 1:length(x$truep)
  }
  # average sample size
  n.average <- weighted.mean(x$ssize, x$prob)
  n.min <- min(x$ssize[x$prob>0])
  n.max <- max(x$ssize[x$prob>0])
  tab0 <- cbind(n.average, n.min, n.max)
  rownames(tab0) <- "Sample size"
  colnames(tab0) <- c("Mean", "Minimum", "Maximum")
  exp <- c(NA, x$exp)
  rec <- xtabs(x$prob~x$mtd)
  tab <- signif(rbind(exp, rec), 3)
  rownames(tab) <- c("Experimentation proportion", "Recommendation proportion")
  colnames(tab) <- c(paste("<", dose[1]), dose)
  names(dimnames(tab)) <- c("", "Doses")
  if(is.null(tox.cutpoints)){
    tox.cutpoints <- seq(0, 1, by=0.2)
  } else {
    tox.cutpoints <- unique(c(0, tox.cutpoints, 1))
  }
  exp.tox <- xtabs(x$exp~cut(x$truep, tox.cutpoints, include.lowest=T))
  rec.tox <- xtabs(rec[-1]~cut(x$truep, tox.cutpoints, include.lowest=T))
  tab2 <- signif(rbind(exp.tox, rec.tox), 3)
  rownames(tab2) <- c("Experimentation proportion", "Recommendation proportion*")
  names(dimnames(tab2)) <- c("", "Probability of DLT")

  tab3 <- signif(rbind(x$n.average, x$dlt.average), 3)
  rownames(tab3) <- c("Average number of patients", "Average number of DLTs")
  colnames(tab3) <- dose
  names(dimnames(tab3)) <- c("", "Doses")

  print(tab0)
  cat("\n")
  print(tab)
  cat("\n")
  print(tab2)
  cat("\n * Amongst those trials that recommend an MTD\n")
  cat("\n")
  print(tab3)
}

#####################################################################
#
# threep3 - Generates all possible 3+3 trial pathways,  probabilities
#            of occurrence,  MTD recommendation and sample size
# AUTHOR: Graham Wheeler
#
# ARGUMENTS
# truep - vector of probabilities for DLT at each dose level
#              to be investigated
# start - starting dose level for the design (defaults to 1)
#    dose    --> optional vector of dose labels (for printing and plotting purposes)
#
# VALUES
#       "prob" - probability of each trial occurring
#       "ssize" - sample size per trial
#       "mtd" - final MTD recommendation per trial
#     "exp" - experimentation proportions across all possible 3+3 trials
#       "truep" - true vector of probabilities
#
####################################################################
#####################################################################
#
# threep3 - Generates all possible 3+3 trial pathways,  probabilities
#            of occurrence,  MTD recommendation and sample size
# AUTHORS: Graham Wheeler,  Michael Sweeting
#
# ARGUMENTS
# truep - vector of probabilities for DLT at each dose level
#              to be investigated
# start - starting dose level for the design (defaults to 1)
#    dose    --> optional vector of dose labels (for printing and plotting purposes)
#
# VALUES
#       "prob" - probability of each trial occurring
#       "ssize" - sample size per trial
#       "mtd" - final MTD recommendation per trial
#     "exp" - experimentation proportions across all possible 3+3 trials
#       "truep" - true vector of probabilities
#
####################################################################



#' Calculate all possible trial pathways for the standard 3+3 design,  together
#' with their probability of occurring
#'
#' All possible pathways of a standard 3+3 design are calculated and assigned a
#' probability of occurring. This facilitates the calculation of operating
#' characteristics,  using \code{\link{print.threep3}} and
#' \code{\link{plot.threep3}}.
#'
#' The dose-escalation schema used here relates to that defined by Chang et al.
#' (2006),  which incorporates dose de-escalation also. Variations of this
#' design exist in the literature (see Storer 1989,  Reiner et al. 1999).
#'
#' The first cohort of three patients are administered the starting dose
#' (usually the lowest dose). The trial then proceeds as follows: \itemize{
#' \item If none of the three patients experience a DLT,  then dose the next
#' three patients at the next highest dose level; \item If one of the three
#' patients last treated experiences a DLT,  then dose the next three patients
#' at the current dose level; \item If at least two patients in the first dose
#' level experience a DLT the trial is stopped for safety and no dose is
#' recommended; } Escalation / de-escalation rules to the next dose level for
#' subsequent cohorts proceed as follows: \itemize{ \item Escalate: If 0/3 or
#' at most 1/6 DLTs are observed in the current cohort AND the next highest
#' dose has not yet been tested; \item Stay at current dose level: If 1/3 DLTs
#' have been observed at this level. Dose a further three patients at the same
#' level; \item De-Escalate: If at least two out of three to six patients
#' experience DLTs at the current dose level AND fewer than six patients have
#' been dosed at the next lowest level }
#'
#' If none of the rules above are satisfied then the trial stops. If the
#' current dose level has at most one DLT observed then this is claimed to be
#' the MTD,  otherwise the dose level below is deemed to be the MTD.
#'
#' If dose-escalation extends to doses outside of that defined by \code{dose},
#' the MTD is determined to be the largest dose in \code{dose}.
#'
#' @param truep A vector of length \code{k} (the number of doses being
#' considered in the trial),  with values equal to the true probabilities of
#' toxicity at the dose levels.
#' @param start Starting dose level. Defaults to 1,  i.e. the lowest dose level
#' @param dose Optional vector of length \code{k} of actual doses for
#' presentation purposes
#' @return \code{threep3} returns an object of class "threep3". The function
#' \code{\link{print}} (i.e. \code{\link{print.threep3}}) can be used to obtain
#' operating characteristics of the design used.
#'
#' An object of class "threep3" is a list with the following components:
#' \item{prob}{A vector with the probabilities of each design occurring. As all
#' possible designs are calculated,  this vector sums to one} \item{ssize}{A
#' vector with the sample size of each design} \item{mtd}{A vector of dose
#' levels giving the recommended maximum tolerated dose (MTD) at the end of the
#' trial} \item{exp}{A vector of length \code{k} giving the average trial
#' experimentation proportions at each dose level} \item{dlt.no}{A vector with
#' the number of toxicities (DLTs) that occur in each trial} \item{truep}{The
#' true probabilities of toxicity at each dose level,  specified by the user}
#' \item{dose}{The actual doses as supplied in the function arguments}
#' \item{n.average}{The average number of patients dosed at each level}
#' \item{dlt.average}{The average number of DLTs experienced at each dose
#' level} \item{all.designs}{A matrix containing all possible 3+3 designs,  with
#' each row representing a different design. Columns labelled "d k" and "tox k"
#' represent the dose level and number of toxicities for the kth cohort,
#' respectively.}
#' @author Graham Wheeler \email{graham.wheeler@@mrc-bsu.cam.ac.uk} (MRC
#' Biostatistics Unit,  Cambridge,  UK) and
#'
#' Michael Sweeting \email{mjs212@@medschl.cam.ac.uk} (University of Cambridge,
#' UK)
#' @seealso \code{\link{threep3}}
#' @references Sweeting M.,  Mander A.,  Sabin T. \pkg{bcrm}: Bayesian Continual
#' Reassessment Method Designs for Phase I Dose-Finding Trials. \emph{Journal
#' of Statistical Software} (2013) 54: 1--26.
#' \url{http://www.jstatsoft.org/article/view/v054i13}
#'
#' Chang A.,  Ganz P.,  Hayes D.,  Kinsella T.,  Pass H.,  Schiller J.,  Stone R.,
#' Strecher V. \emph{Oncology: An Evidence-Based Approach}. Springer (2006).
#'
#' Storer B. Design and Analysis of Phase I Clinical Trials. \emph{Biometrics}
#' (1989) 45: 925--937.
#'
#' Reiner E.,  Paoletti X.,  O'Quigley J. Operating characteristics of the
#' standard phase I clinical trial design. \emph{Computational Statistics &
#' Data Analysis} (1999) 30: 303--315.
#'
#' Neuenschwander B.,  Branson M.,  Gsponer T. Critical aspects of the Bayesian
#' approach to phase I cancer trials. \emph{Statistics in Medicine} (2008) 27:
#' 2420--2439.
#' @examples
#'
#' ## What are the operating characteristics of a standard 3+3 design if we conside only the first
#' ## 12 doses of the dose-escalation cancer trial example as described in Neuenschwander et al 2008.
#' ## Pre-defined doses
#' dose <- c(1, 2.5, 5, 10, 15, 20, 25, 30, 40, 50, 75, 100)
#' ## Pre-specified probabilities of toxicity
#' p.tox0 <- c(0.010, 0.015, 0.020, 0.025, 0.030, 0.040, 0.050, 0.100, 0.170, 0.300, 0.400, 0.500)
#'
#' \dontrun{
#' design.threep3 <- threep3(p.tox0, dose)
#' print(design.threep3)
#' plot(design.threep3)
#' }
#'
#' @export threep3
threep3 <- function(truep, start=1, dose=NULL){
  # Check that dose is the same length as truep
  if(!is.null(dose) & length(dose)!=length(truep)) stop("Length of 'dose' must be the same as the length of 'truep'.")

  # Define number of doses and max. cohort number
  # 'mcplus1' used in computation
  doses <- length(truep)
  mcohort <- 2*doses
  mcplus1 <- mcohort+1

  # Check that starting dose is not top dose
  if(start==doses) stop("Starting dose cannot be uppermost dose")

  # Begin deriving pathways
  pmat <- as.data.frame(matrix(NA, nrow=1, ncol=2*mcplus1+1))
  colnames(pmat) <- c("stop", "desc", paste(c("d", "tox"), rep(1:mcohort, each=2)), paste("d", mcplus1))
  pmat[1, 1:3] <- c(0, 0, start)
  pmat <- pmat[rep(seq_len(nrow(pmat)),  rep(4, nrow(pmat))), ]
  pmat[, "tox 1"] <- c(0, 1, 2, 3)

  pmat[pmat[, "tox 1"]==0, "d 2"] <- start+1
  pmat[pmat[, "tox 1"]==1, "d 2"] <- start
  pmat[pmat[, "tox 1"]>1 & start>1, "d 2"] <- start-1
  pmat[pmat[, "tox 1"]>1 & start==1, "stop"] <- 1
  pmat[pmat[, "tox 1"]>1, "desc"] <- 1
  stopped.pmat <- pmat[pmat$stop==1, -2]
  all.designs <- stopped.pmat
  prob <- ssize <- mtd <- dlt.no <- NULL
  exp <- 0

  if(start==1){
    # Probabilties of stopped 3+3 trials occuring
    dose.mat <- stopped.pmat[, grep("d", names(stopped.pmat))]
    tox.mat <- stopped.pmat[, grep("tox", names(stopped.pmat))]
    prob <- apply(matrix(dbinom(as.matrix(tox.mat), 3, truep[as.matrix(dose.mat)]), nrow=nrow(dose.mat)), 1, prod, na.rm=T)

    # Calculate sample size
    ssize <- 3*apply(!is.na(stopped.pmat[, grep("d", names(stopped.pmat))]), 1, sum)

    # Determine MTD per stopped trial
    last.cohort <- apply(!is.na(stopped.pmat[, grep("d", names(stopped.pmat))]), 1, sum)
    last.drug.column <- paste("d", last.cohort)
    last.drug <- sapply(1:nrow(stopped.pmat), function(j){stopped.pmat[j, last.drug.column[j]]})
    previous.drug.column <- paste("d", last.cohort-1)
    previous.drug <- sapply(1:nrow(stopped.pmat), function(j){ifelse(previous.drug.column[j]=="d 0", 0, stopped.pmat[j, previous.drug.column[j]])})
    last.tox.column <- paste("tox", last.cohort)
    last.tox <- sapply(1:nrow(stopped.pmat), function(j){stopped.pmat[j, last.tox.column[j]]})
    mtd <- rep(NA, nrow(stopped.pmat))
    mtd[last.tox==0] <- last.drug[last.tox==0]
    mtd[last.tox==1 & previous.drug==last.drug] <- last.drug[last.tox==1 & previous.drug==last.drug]-1
    mtd[last.tox==1 & previous.drug!=last.drug] <- last.drug[last.tox==1 & previous.drug!=last.drug]
    mtd[last.tox>1] <- last.drug[last.tox>1]-1

    # Prob. that each dose is experimented on and trial occurs
    exp <- sapply(1:doses, function(j){sum(3*(stopped.pmat[, grep("d", names(stopped.pmat))]==j)*prob/ssize, na.rm=T)})
    # Average number of people dosed at each level
    n.average <- sapply(1:doses, function(j){sum(3*(stopped.pmat[, grep("d", names(stopped.pmat))]==j)*prob, na.rm=T)})

    # Number of subjects who have DLT per trial
    dlt.no <- apply(stopped.pmat[, grep("tox", names(stopped.pmat))], 1, sum, na.rm=T)
    # Average number of DLTs at each dose
    dlt.average <- sapply(1:doses, function(j){sum((stopped.pmat[, grep("tox", names(stopped.pmat))]*prob)[stopped.pmat[, grep("d", names(stopped.pmat))]==j], na.rm=T)})


  }

  for(i in 3:mcplus1){
    cat(paste(round(100*i/mcplus1), "% complete\n", sep=""))
    dd <- as.character(paste("d", i))
    td <- as.character(paste("tox", i))
    dc <- as.character(paste("d", i-1))
    tc <- as.character(paste("tox", i-1))
    db <- as.character(paste("d", i-2))
    tb <- as.character(paste("tox", i-2))

    ## Creates new data.frame with 1,  2,  or 3 toxicities for every continued trial
    pmat <- pmat[rep(which(pmat[, "stop"]==0),  each=4), ]
    pmat[, tc] <- 0:3

    pmat[pmat[, tc]==0 & pmat[, "desc"]==0 & pmat[, dc]+1 <=doses, dd] <-  pmat[pmat[, tc]==0 & pmat[, "desc"]==0 & pmat[, dc]+1 <=doses, dc]+1
    pmat[pmat[, tc]==1 & pmat[, "desc"]==0 & pmat[, tb]==0, dd] <-  pmat[pmat[, tc]==1 & pmat[, "desc"]==0 & pmat[, tb]==0, dc]
    pmat[pmat[, tc]==1 & pmat[, "desc"]==0 & pmat[, tb]==1 & pmat[, dc]-1 >=1 , dd] <-  pmat[pmat[, tc]==1 & pmat[, "desc"]==0 & pmat[, tb]==1, dc]-1
    pmat[pmat[, tc]==0 & pmat[, "desc"]==1 , dd] <-  pmat[pmat[, tc]==0 & pmat[, "desc"]==1, dc]
    pmat[pmat[, tc]==1 & pmat[, "desc"]==1,  dd] <-  pmat[pmat[, tc]==1 & pmat[, "desc"]==1, dc]
    pmat[pmat[, tc]>1 & pmat[, dc]-1 >=1, "desc"] <-  1

    pmat[pmat[, tc]==1 & pmat[, "desc"]==0 & pmat[, tb]==1 & pmat[, dc]-1 >=1, "desc"] <-  1
    pmat[pmat[, tc]>1 & pmat[, dc]-1 >=1, dd] <-  pmat[pmat[, tc]>1 & pmat[, dc]-1>=1 , dc]-1


    excluding.dd <- names(pmat)[grepl("d ", names(pmat)) & names(pmat)!=dd]
    cnt <- apply(pmat[!is.na(pmat[, dd]), dd]==pmat[!is.na(pmat[, dd]), excluding.dd], 1, sum, na.rm=T)

    pmat[!is.na(pmat[, dd]), dd][cnt>1] <- NA
    pmat[is.na(pmat[, dd]), "stop"] <- 1
    stopped.pmat <- pmat[pmat$stop==1, -2]
    all.designs <- rbind(all.designs, stopped.pmat)
    if(dim(stopped.pmat)[1]>0){
      # Probabilties of stopped 3+3 trials occuring
      dose.mat <- stopped.pmat[, grep("d", names(stopped.pmat))[1:(i-1)]]
      tox.mat <- stopped.pmat[, grep("tox", names(stopped.pmat))[1:(i-1)]]


      # Add these probabilities to prob
      prob.new <- apply(matrix(dbinom(as.matrix(tox.mat), 3, truep[as.matrix(dose.mat)]), nrow=nrow(dose.mat)), 1, prod, na.rm=T)
      prob <- c(prob, prob.new)

      # Calculate sample size and determine MTD per stopped trial
      # Add them to existing ssize and mtd vectors
      ssize.new <- rep(3*(i-1), nrow(stopped.pmat))
      ssize <- c(ssize, ssize.new)
      last.drug <- stopped.pmat[, dc]
      previous.drug <- stopped.pmat[, db]
      last.tox <- stopped.pmat[, tc]

      mtd.new <- rep(NA, nrow(stopped.pmat))
      mtd.new[last.tox==0] <- last.drug[last.tox==0]
      mtd.new[last.tox==1 & previous.drug==last.drug] <- last.drug[last.tox==1 & previous.drug==last.drug]-1
      mtd.new[last.tox==1 & previous.drug!=last.drug] <- last.drug[last.tox==1 & previous.drug!=last.drug]
      mtd.new[last.tox>1] <- last.drug[last.tox>1]-1
      mtd <- c(mtd, mtd.new)

      # Prob. that each dose is experimented on and trial occurs (summing over all trials)
      exp <- exp+sapply(1:doses, function(j){sum(3*(stopped.pmat[, grep("d", names(stopped.pmat))]==j)*prob.new/ssize.new, na.rm=T)})
      # Average number of people dosed at each level
      n.average <- n.average+sapply(1:doses, function(j){sum(3*(stopped.pmat[, grep("d", names(stopped.pmat))]==j)*prob.new, na.rm=T)})


      # Number of subjects who have DLT per trial
      dlt.no <- c(dlt.no, apply(stopped.pmat[, grep("tox", names(stopped.pmat))], 1, sum, na.rm=T))
      # Average number of DLTs at each dose
      dlt.average <- dlt.average+sapply(1:doses, function(j){sum((stopped.pmat[, grep("tox", names(stopped.pmat))]*prob.new)[stopped.pmat[, grep("d", names(stopped.pmat))]==j], na.rm=T)})


    }
  }
  obj <- list(prob=prob, ssize=ssize, mtd=mtd, exp=exp, dlt.no=dlt.no, truep=truep, dose=dose, n.average=n.average, dlt.average=dlt.average, all.designs=all.designs)
  class(obj) <- "threep3"
  return(obj)
}

