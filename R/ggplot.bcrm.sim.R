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
plot.bcrm.sim <- function(x, trajectories=FALSE, file=NULL, threep3=FALSE, ...){
  dose <- if(is.null(x[[1]]$dose)) x[[1]]$sdose else x[[1]]$dose
  dose.label <- if(is.null(x[[1]]$dose)) "Standardised dose" else "Dose"

  if(trajectories){
    ## sample size
    n <- sapply(x, function(i){dim(i$data)[1]})
    traj.mat <- matrix(NA, nrow=length(x), ncol=max(n))
    tox.mat <- matrix(NA, nrow=length(x), ncol=max(n))
    for(i in 1:length(x)){
      traj.mat[i, 1:n[i]] <- x[[i]]$data$dose
      tox.mat[i, 1:n[i]] <- x[[i]]$data$tox
    }
    traj.df <- data.frame(patient=rep(1:max(n), each=5), Statistic=factor(rep(c("Minimum", "Lower Quartile", "Median", "Upper Quartile", "Maximum"), max(n)), levels=c("Minimum", "Lower Quartile", "Median", "Upper Quartile", "Maximum")), traj=c(apply(traj.mat, 2, quantile, na.rm=T)))
    df <- data.frame(patient=rep(1:max(n), each=length(x)), sim=rep(1:length(n), max(n)), traj=c(traj.mat), Toxicity=ifelse(c(tox.mat)==1, "Yes", "No"))
    lt <- c("Median" = 1, "Lower Quartile" = 2, "Upper Quartile" = 2,  "Minimum" = 4, "Maximum"=4)
    cols <-  c("No" = "black", "Yes" = "red")
    if(length(x) > 1){
      #a <- ggplot()+geom_point(aes(x=patient, y=traj, col=Toxicity), data=df, position="jitter", alpha=0.2)+
      #  scale_colour_manual(values=cols)+
      #  xlab("Patient")+ylab("Dose Level")+ggtitle("Experimentation and Toxicities")
      b <- ggplot() +
        geom_step(aes(x=patient, y=traj, group=Statistic, linetype=Statistic), size=1.2, colour="blue", data=traj.df) +
        scale_linetype_manual(values=lt) +
        xlab("Patient")+ylab("Dose Level")
      if(!is.null(file))
        pdf(paste(file, ".pdf", sep=""), ...)
      #grid.newpage()
      #pushViewport(viewport(layout=grid.layout(2, 1)))
      #vplayout <- function(x, y)  viewport(layout.pos.row=x, layout.pos.col=y)
      #print(a, vp=vplayout(1, 1))
      print(b)
      if(!is.null(file))
        dev.off()
    } else {
      a <- ggplot()+scale_colour_manual(values=cols)+
        xlab("Patient")+ylab("Dose Level")+
        geom_point(aes(x=patient, y=traj, col=Toxicity), data=df)
      print(a)
      if(!is.null(file))
        ggsave(paste(file, ".pdf", sep=""), ...)
    }
  } else {
    if(threep3 & is.null(x[[1]]$threep3)){
      cat("\n Calculating 3+3 operating characteristics....\n")
      x[[1]]$threep3 <- threep3(x[[1]]$truep, x[[1]]$start)
    }
    # sample size
    n <- sapply(x, function(i){dim(i$data)[1]})
    df.n <- data.frame(n)
    if(!threep3){
      a <- if(min(n)!=max(n)){
        ggplot()+geom_histogram(aes(x=n, y=100*..density..), data=df.n, binwidth=1)+xlab("Sample size")+ylab("Percent")+ggtitle("Sample size")
      } else { NULL }
    } else {
      n.threep3 <- x[[1]]$threep3$ssize
      df.n.threep3 <- data.frame(n=c(n, n.threep3),
                                 weight=c(rep(1/length(n), length(n)), data[[1]]$threep$prob),
                                 Method=rep(c("CRM", "3+3"), c(length(n), length(n.threep3))))
      a <- ggplot() +
        stat_bin(aes(x=n, y=100*..density.., weight=weight, fill=Method), data=df.n.threep3, binwidth=1, position="dodge") +
        xlab("Sample size") +
        ylab("Percent") +
        ggtitle("Sample size")
    }

    # experimentation
    exp <- rep(dose, apply(sapply(x, function(i){ (i$tox + i$notox)} ), 1, sum))
    if(!threep3){
      df.exp <- data.frame(exp=factor(exp))
      b <- ggplot() +
        geom_bar(aes(x=exp, y=100*..count../sum(..count..)), data=df.exp) +
        xlab(dose.label) +
        ylab("Percent") +
        ggtitle("Experimentation")
    } else {
      exp.threep3 <- rep(dose, 10000*x[[1]]$threep3$exp)
      df.exp.threep3 <- data.frame(exp=factor(c(exp, exp.threep3)),
                                   Method=rep(c("CRM", "3+3"), c(length(exp), length(exp.threep3))),
                                   weight=c(rep(1/length(exp), length(exp)), rep(1/length(exp.threep3), length(exp.threep3))))
      b <- ggplot() +
        geom_bar(aes(x=exp, y=100*..count.., weight=weight, fill=Method), data=df.exp.threep3, position="dodge") +
        xlab(dose.label) +
        ylab("Percent") +
        ggtitle("Experimentation")
    }

    # recommendation
    rec <- dose[sapply(x, function(i){i$ndose[[length(i$ndose)]]$ndose})]
    if(!threep3){
      df.rec <- data.frame(rec=factor(rec))
      c <- ggplot() +
        geom_bar(aes(x=rec, y=100*..count../sum(count)), data=df.rec) +
        xlab(dose.label) +
        ylab("Percent") +
        ggtitle("Recommendation")
    } else {
      rec.threep3 <- dose[x[[1]]$threep3$mtd]
      df.rec.threep3 <- data.frame(rec=factor(c(rec, rec.threep3)), weight=c(rep(1/length(rec), length(rec)), x[[1]]$threep$prob[x[[1]]$threep3$mtd!=0]), Method=rep(c("CRM", "3+3"), c(length(rec), length(rec.threep3))))
      c <- ggplot() +
        geom_bar(aes(x=rec, y=100*..count.., weight=weight, fill=Method), data=df.rec.threep3, position="dodge") +
        xlab(dose.label) +
        ylab("Percent") +
        ggtitle("Recommendation")
    }

    # observed DLTs
    obs <- sapply(x, function(i){100*sum(i$tox)/sum(i$tox+i$notox)})
    if(!threep3){
      bw <- max(diff(range(obs))/30, 1)
      df.obs <- data.frame(obs=obs, bw=bw)
      d <- ggplot()+geom_histogram(aes(bw=bw, x=obs, y=100*..density..*bw), data=df.obs, binwidth=bw)+xlab("Percentage of subjects with DLTs")+ylab("Percent")+ggtitle("DLTs")
    } else {
      obs.threep3 <- 100*x[[1]]$threep3$dlt.no/x[[1]]$threep3$ssize
      bw <- diff(range(c(obs, obs.threep3)))/30
      df.obs.threep3 <- data.frame(bw=bw, obs=c(obs, obs.threep3), weight=c(rep(1/length(obs), length(obs)), x[[1]]$threep$prob), Method=rep(c("CRM", "3+3"), c(length(obs), length(obs.threep3))))
      df.obs.threep3 <- subset(df.obs.threep3, df.obs.threep3$weight>0)
      d <- ggplot()+geom_histogram(aes(bw=bw, x=obs, y=100*..density..*bw, weight=weight, fill=Method), data=df.obs.threep3, binwidth=bw, position="dodge")+xlab("Percentage of subjects with DLTs")+ylab("Percent")+ggtitle("DLTs")
    }
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
}
