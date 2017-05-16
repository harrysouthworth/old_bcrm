#' Helper for getting dose info from a bcrm object
bcrm.get.dose <- function(x){
  dose <- if(is.null(x[[1]]$dose)){
    x[[1]]$sdose
    } else{
      x[[1]]$dose
    }
  attr(dose, "dose.label") <- if(is.null(x[[1]]$dose)) "Standardized dose" else "Dose"

  dose
}

# Next function from https://github.com/tidyverse/ggplot2/wiki/share-a-legend-between-two-ggplot2-graphs
grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)

  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = grid:::unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = grid:::unit.c(unit(1, "npc") - lwidth, lwidth)))

  grid:::grid.newpage()
  grid:::grid.draw(combined)

  # return gtable invisibly
  invisible(combined)

}

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
#' @param data An object of class "bcrm.sim",  as returned by \code{\link{bcrm}}
#' when conducting a simulation.
#' @param mapping Not used.
#' @param ncol, nrow The numbers of rows and columns. Depending on other arguments,
#'   either 3 or 4 plots will be produced, so \code{ncol * nrow} needs to be either
#'   3 or 4.
#' @param trajectories Should a summary plot of the trajectories of
#' administered dose levels be plotted? Defaults to FALSE.
#' @param threep3 Should operating characteristics of a standard 3+3 rule-based
#'   design be plotted alongside the \code{bcrm} design? Defaults to
#'   \code{FALSE}.
#' @param legend.position Position for the legend to go. Should be either
#'   \code{legend.position = "right"} or \code{legend.position = "bottom"} (the
#'   default).
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
ggplot.bcrm.sim <- function(data, mapping=aes(), ncol=2, nrow=2, trajectories=FALSE, threep3=FALSE, quietly=TRUE,
                            legend.position="bottom", ...){
  if(threep3 & is.null(data[[1]]$threep3)){
    if (!quietly){
      message("Calculating 3+3 operating characteristics....")
    }
    data[[1]]$threep3 <- threep3(data[[1]]$truep, data[[1]]$start, quietly=quietly)
  }


  if (trajectories){
    ggplot.bcrm.sim_traj(data, mapping=NULL, threep3=threep3)
  } else {
    a <- ggplot.bcrm.sim_threep3(data, mapping=NULL, threep3=threep3, quietly=quietly)

    b <- ggplot.bcrm.sim_experimentation(data, mapping=NULL, threep3=threep3) +
      theme(legend.position=legend.position)

    cc <- ggplot.bcrm.sim_recommendation(data, mapping=NULL, threep3=threep3)

    d <- ggplot.bcrm.sim_obsDLT(data, mapping=NULL, threep3=threep3)

    if (is.null(a)){
      grid.arrange(b, cc, d, ncol=ncol)
    } else {
      grid_arrange_shared_legend(a, b, cc, d, ncol=ncol, nrow=nrow, position=legend.position)
    }
  }
}


ggplot.bcrm.sim_traj <- function(data, mapping=aes(), threep3=FALSE, ...){
  dose <- bcrm.get.dose(data)

  ## sample size
  n <- sapply(data, function(i){ dim(i$data)[1] })
  traj.mat <- matrix(NA, nrow=length(data), ncol=max(n))
  tox.mat <- matrix(NA, nrow=length(data), ncol=max(n))
  for(i in 1:length(data)){
    traj.mat[i, 1:n[i]] <- data[[i]]$data$dose
    tox.mat[i, 1:n[i]] <- data[[i]]$data$tox
  }
  traj.df <- data.frame(patient=rep(1:max(n), each=5),
                        Statistic=factor(rep(c("Minimum", "Lower Quartile", "Median", "Upper Quartile", "Maximum"), max(n)),
                                         levels=c("Minimum", "Lower Quartile", "Median", "Upper Quartile", "Maximum")),
                        traj=c(apply(traj.mat, 2, quantile, na.rm=TRUE)))
  df <- data.frame(patient = rep(1:max(n), each=length(data)),
                   sim = rep(1:length(n), max(n)), traj=c(traj.mat),
                   Toxicity = ifelse(c(tox.mat)==1, "Yes", "No"))

  lt <- c("Median" = 1, "Lower Quartile" = 2, "Upper Quartile" = 2,  "Minimum" = 4, "Maximum"=4)
  cols <-  c("No" = "black", "Yes" = "red")

  if(length(data) > 1){
    ggplot(data=traj.df, aes(x=patient, y=traj, group=Statistic, linetype=Statistic)) +
      geom_step(size=1.2, colour="blue") +
      scale_linetype_manual(values=lt) +
      xlab("Patient") +
      ylab("Dose Level")
  } else {
    ggplot(data=df, aes(x=patient, y=traj, col=Toxicity)) +
      scale_colour_manual(values=cols) +
      xlab("Patient") +
      ylab("Dose Level") +
      geom_point()
  }
}

ggplot.bcrm.sim_threep3 <- function(data, mapping=aes(), threep3=FALSE, ...){
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


ggplot.bcrm.sim_experimentation <- function(data, mapping=aes(), threep3=FALSE, ...){
  dose <- bcrm.get.dose(data)

  # experimentation
  exp <- rep(dose, apply(sapply(data, function(i){ (i$tox + i$notox)} ), 1, sum))

  if(!threep3){
    df.exp <- data.frame(exp=factor(exp))

    ggplot(data=df.exp, aes(x=exp, y = 100 * ..count.. / sum(..count..))) +
      geom_bar() +
      xlab(attributes(dose)[["dose.label"]]) +
      ylab("Percent") +
      ggtitle("Experimentation")
  } else {
    exp.threep3 <- rep(dose, 10000 * data[[1]]$threep3$exp)

    df.exp.threep3 <- data.frame(exp=factor(c(exp, exp.threep3)),
                                 Method=rep(c("CRM", "3+3"), c(length(exp), length(exp.threep3))),
                                 weight=c(rep(1/length(exp), length(exp)), rep(1/length(exp.threep3), length(exp.threep3))))

    ggplot(data=df.exp.threep3, aes(x=exp, y=100*..count.., weight=weight, fill=Method)) +
      geom_bar(position="dodge") +
      xlab(attributes(dose)[["dose.label"]]) +
      ylab("Percent") +
      ggtitle("Experimentation")
  }
}

ggplot.bcrm.sim_recommendation <- function(data, mapping=aes(), threep3=FALSE, ...){
  dose <- bcrm.get.dose(data)

  # recommendation
  rec <- dose[sapply(data, function(i){ i$ndose[[length(i$ndose)]]$ndose })]

  if(!threep3){
    df.rec <- data.frame(rec=factor(rec))
    ggplot(data=df.rec, aes(x=rec, y = 100 * ..count.. / sum(count))) +
      geom_bar() +
      xlab(attributes(dose)[["dose.label"]]) +
      ylab("Percent") +
      ggtitle("Recommendation")
  } else {
    rec.threep3 <- dose[data[[1]]$threep3$mtd]
    df.rec.threep3 <- data.frame(rec=factor(c(rec, rec.threep3)),
                                 weight=c(rep(1/length(rec), length(rec)), data[[1]]$threep$prob[data[[1]]$threep3$mtd!=0]),
                                 Method=rep(c("CRM", "3+3"), c(length(rec), length(rec.threep3))))
    ggplot(data=df.rec.threep3, aes(x=rec, y=100*..count.., weight=weight, fill=Method)) +
      geom_bar(position="dodge") +
      xlab(attributes(dose)[["dose.label"]]) +
      ylab("Percent") +
      ggtitle("Recommendation")
  }
}

ggplot.bcrm.sim_obsDLT <- function(data, mapping=aes(), threep3=FALSE, ...){
  dose <- bcrm.get.dose(data)

  # observed DLTs
  obs <- sapply(data, function(i){ 100 * sum(i$tox) / sum(i$tox + i$notox)})
  if(!threep3){
    bw <- max(diff(range(obs)) / 30, 1)
    df.obs <- data.frame(obs=obs, bw=bw)

    ggplot(data=df.obs, aes(bw=bw, x=obs, y=100*..density..*bw)) +
      geom_histogram(binwidth=bw) +
      xlab("Percentage of subjects with DLTs") +
      ylab("Percent") +
      ggtitle("DLTs")
  } else {
    obs.threep3 <- 100 * data[[1]]$threep3$dlt.no / data[[1]]$threep3$ssize
    bw <- diff(range(c(obs, obs.threep3))) / 30
    df.obs.threep3 <- data.frame(bw=bw,
                                 obs=c(obs, obs.threep3),
                                 weight=c(rep(1/length(obs), length(obs)), data[[1]]$threep$prob),
                                 Method=rep(c("CRM", "3+3"), c(length(obs), length(obs.threep3))))

    df.obs.threep3 <- subset(df.obs.threep3, df.obs.threep3$weight > 0)

    ggplot(data=df.obs.threep3, aes(bw=bw, x=obs, y=100*..density..*bw, weight=weight, fill=Method)) +
      geom_histogram(binwidth=bw, position="dodge") +
      xlab("Percentage of subjects with DLTs") +
      ylab("Percent") +
      ggtitle("DLTs")
  }
}
