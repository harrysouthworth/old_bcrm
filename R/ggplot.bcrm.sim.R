#-----------------------------------------------------------------------
#    Plot function for an object of class bcrm.sim
#    threep3     --> If TRUE (default is FALSE) then operating characteristics of the simulated design are compared against a standard rule-based 3+3 design
# -----------------------------------


#' Plot the operating characteristics from the simulated trials
#'
#' Plots of the operating characteristics obtained from a CRM simulation.
#'
#' This function plots the sample size distribution (if variable),  the
#' experimentation distribution,  the recommended dose distribution and the
#' percentage of subjects who experience the toxicity outcome (dose-limiting
#' toxicity). If \code{trajectories = TRUE} then summary statistics of
#' administered dose levels for each patient are plotted instead. If
#' \code{threep3 = TRUE} then the operating characteristics of the standard 3+3
#' design are plotted alongside those of the \code{bcrm} design (see
#' \code{\link{threep3}} for more details).
#'
#' @param x An object of class "bcrm.sim",  as returned by \code{\link{bcrm}}
#' when conducting a simulation.
#' @param trajectories Should a summary plot of the trajectories of
#' administered dose levels be plotted? Defaults to FALSE.
#' @param file File name where the operating characteristic plot is stored,  in
#' a pdf format.
#' @param threep3 Should operating characteristics of a standard 3+3 rule-based
#' design be plotted alongside the \code{bcrm} design? Defaults to
#' \code{FALSE}.
#' @param ... Further arguments passed to or from other methods
#' @author Michael Sweeting \email{mjs212@@medschl.cam.ac.uk} (University of
#' Cambridge,  UK)
#' @seealso \code{\link{print.bcrm.sim}},  \code{\link{bcrm}},
#' \code{\link{threep3}}
#' @references Sweeting M.,  Mander A.,  Sabin T. \pkg{bcrm}: Bayesian Continual
#' Reassessment Method Designs for Phase I Dose-Finding Trials. \emph{Journal
#' of Statistical Software} (2013) 54: 1--26.
#' \url{http://www.jstatsoft.org/article/view/v054i13}
#' @method ggplot bcrm.sim
#' @export
ggplot.bcrm.sim <- function(data, mapping=aes(), trajectories=FALSE, file=NULL, threep3=FALSE, quietly=TRUE, ...){
  dose <- if(is.null(data[[1]]$dose)) data[[1]]$sdose else data[[1]]$dose
  dose.label <- if(is.null(data[[1]]$dose)) "Standardized dose" else "Dose"

  if (trajectories){
    ggplot.bcrm.sim_traj(data, mapping=NULL)
  } else {
    a <- ggplot.bcrm.sim_threep3(data, mapping=NULL)

    b <- ggplot.bcrm.sim_experimentation(data, mapping=NULL)

    cc <- ggplot.bcrm.sim_recommendation(data, mapping=NULL)

    d <- ggplot.bcrm.sim_obsDLT(data, mapping=NULL)


    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout=grid::grid.layout(2, 2)))
    vplayout <- function(x, y)  grid::viewport(layout.pos.row=x, layout.pos.col=y)
    if(!is.null(a)) print(a, vp=vplayout(1, 1))
    print(b, vp=vplayout(1, 2))
    print(c, vp=vplayout(2, 1))
    print(d, vp=vplayout(2, 2))
  }
}


ggplot.bcrm.sim_traj <- function(data, mapping=aes(), ...){
  ## sample size
  n <- sapply(data, function(i){ dim(i$data)[1] })
  traj.mat <- matrix(NA, nrow=length(data), ncol=max(n))
  tox.mat <- matrix(NA, nrow=length(data), ncol=max(n))
  for(i in 1:length(data)){
    traj.mat[i, 1:n[i]] <- x[[i]]$data$dose
    tox.mat[i, 1:n[i]] <- x[[i]]$data$tox
  }
  traj.df <- data.frame(patient=rep(1:max(n), each=5),
                        Statistic=factor(rep(c("Minimum", "Lower Quartile", "Median", "Upper Quartile", "Maximum"), max(n)),
                                         levels=c("Minimum", "Lower Quartile", "Median", "Upper Quartile", "Maximum")),
                        traj=c(apply(traj.mat, 2, quantile, na.rm=T)))
  df <- data.frame(patient=rep(1:max(n), each=length(data)),
                   sim=rep(1:length(n), max(n)), traj=c(traj.mat),
                   Toxicity=ifelse(c(tox.mat)==1, "Yes", "No"))

  lt <- c("Median" = 1, "Lower Quartile" = 2, "Upper Quartile" = 2,  "Minimum" = 4, "Maximum"=4)
  cols <-  c("No" = "black", "Yes" = "red")

  if(length(data) > 1){
    b <- ggplot(data=traj.df, aes(x=patient, y=traj, group=Statistic, linetype=Statistic)) +
      geom_step(size=1.2, colour="blue") +
      scale_linetype_manual(values=lt) +
      xlab("Patient") +
      ylab("Dose Level")

    b
  } else {
    a <- ggplot(data=df, aes(x=patient, y=traj, col=Toxicity)) +
      scale_colour_manual(values=cols) +
      xlab("Patient")+ylab("Dose Level") +
      geom_point()

    a
  }
}

ggplot.bcrm.sim_threep3 <- function(data, mapping=aes(), ...){
  if(threep3 & is.null(data[[1]]$threep3)){
    if (!quietly){
      message("Calculating 3+3 operating characteristics....")
    }
    data[[1]]$threep3 <- threep3(data[[1]]$truep, data[[1]]$start, quietly=quietly)
  }
  # sample size
  n <- sapply(data, function(i){ dim(i$data)[1]} )
  df.n <- data.frame(n)
  if(!threep3){
    if(min(n) != max(n)){
      ggplot(data=df.n, aes(x=n, y=100*..density..)) +
        geom_histogram(binwidth=1) +
        xlab("Sample size") +
        ylab("Percent") +
        ggtitle("Sample size")
    } else {
      NULL
    }
  } else { # threep3
    n.threep3 <- data[[1]]$threep3$ssize
    df.n.threep3 <- data.frame(n=c(n, n.threep3),
                               weight=c(rep(1/length(n), length(n)), data[[1]]$threep$prob),
                               Method=rep(c("CRM", "3+3"), c(length(n), length(n.threep3))))
    ggplot(data=df.n.threep3, aes(x=n, y=100*..density.., weight=weight, fill=Method)) +
      stat_bin(binwidth=1, position="dodge") +
      xlab("Sample size") +
      ylab("Percent") +
      ggtitle("Sample size")
  }
}


ggplot.bcrm.sim_experimentation <- function(data, mapping=aes(), ...){
  # experimentation
  exp <- rep(dose, apply(sapply(data, function(i){ (i$tox + i$notox)} ), 1, sum))

  if(!threep3){
    df.exp <- data.frame(exp=factor(exp))

    ggplot(data=df.exp, aes(x=exp, y = 100 * ..count.. / sum(..count..))) +
      geom_bar() +
      xlab(dose.label) +
      ylab("Percent") +
      ggtitle("Experimentation")
  } else {
    exp.threep3 <- rep(dose, 10000 * data[[1]]$threep3$exp)

    df.exp.threep3 <- data.frame(exp=factor(c(exp, exp.threep3)),
                                 Method=rep(c("CRM", "3+3"), c(length(exp), length(exp.threep3))),
                                 weight=c(rep(1/length(exp), length(exp)), rep(1/length(exp.threep3), length(exp.threep3))))

    ggplot(data=df.exp.threep3, aes(x=exp, y=100*..count.., weight=weight, fill=Method)) +
      geom_bar(position="dodge") +
      xlab(dose.label) +
      ylab("Percent") +
      ggtitle("Experimentation")
  }
}

ggplot.bcrm.sim_recommendation(data, mapping=aes(), ...){
  # recommendation
  rec <- dose[sapply(data, function(i){ i$ndose[[length(i$ndose)]]$ndose })]

  if(!threep3){
    df.rec <- data.frame(rec=factor(rec))
    c <- ggplot(data=df.rec, aes(x=rec, y = 100 * ..count.. / sum(count))) +
      geom_bar() +
      xlab(dose.label) +
      ylab("Percent") +
      ggtitle("Recommendation")
  } else {
    rec.threep3 <- dose[data[[1]]$threep3$mtd]
    df.rec.threep3 <- data.frame(rec=factor(c(rec, rec.threep3)),
                                 weight=c(rep(1/length(rec), length(rec)), x[[1]]$threep$prob[x[[1]]$threep3$mtd!=0]),
                                 Method=rep(c("CRM", "3+3"), c(length(rec), length(rec.threep3))))
    c <- ggplot(data=df.rec.threep3, aes(x=rec, y=100*..count.., weight=weight, fill=Method)) +
      geom_bar(position="dodge") +
      xlab(dose.label) +
      ylab("Percent") +
      ggtitle("Recommendation")
  }
}

ggplot.bcrm.sim_obsDLT <- function(data, mapping=aes(), ...){
  # observed DLTs
  obs <- sapply(x, function(i){100*sum(i$tox)/sum(i$tox+i$notox)})
  if(!threep3){
    bw <- max(diff(range(obs))/30, 1)
    df.obs <- data.frame(obs=obs, bw=bw)
    d <- ggplot() +
      geom_histogram(aes(bw=bw, x=obs, y=100*..density..*bw), data=df.obs, binwidth=bw) +
      xlab("Percentage of subjects with DLTs")+ylab("Percent")+ggtitle("DLTs")
  } else {
    obs.threep3 <- 100*x[[1]]$threep3$dlt.no/x[[1]]$threep3$ssize
    bw <- diff(range(c(obs, obs.threep3)))/30
    df.obs.threep3 <- data.frame(bw=bw, obs=c(obs, obs.threep3), weight=c(rep(1/length(obs), length(obs)), x[[1]]$threep$prob), Method=rep(c("CRM", "3+3"), c(length(obs), length(obs.threep3))))
    df.obs.threep3 <- subset(df.obs.threep3, df.obs.threep3$weight>0)
    d <- ggplot(data=df.obs.threep3) +
      geom_histogram(aes(bw=bw, x=obs, y=100*..density..*bw, weight=weight, fill=Method),
                     binwidth=bw, position="dodge") +
      xlab("Percentage of subjects with DLTs") +
      ylab("Percent") +
      ggtitle("DLTs")
  }
}
