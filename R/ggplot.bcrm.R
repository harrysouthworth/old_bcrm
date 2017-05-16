#' Plot the estimated dose-toxicity curve
#'
#' The estimated dose-toxicity curve using the Bayesian continuous reassessment
#' method is plotted for the patients thus far recruited into the trial
#'
#' The estimated 2.5\%,  25\%,  50\%,  75\%,  97.5\% quantiles of the probability
#' of toxicity are plotted for each dose. Additionally,  a histogram of the
#' number of toxicities and non-toxicities is plotted at each experimented
#' dose.
#'
#' If \code{trajectory = TRUE} then the sequential dose trajectory and observed
#' toxicities are plotted.
#'
#' @param x An object of class "bcrm",  as returned by \code{\link{bcrm}}
#' @param each Should posterior summaries be plotted after each recruited
#' cohort? Defaults to FALSE.
#' @param trajectory Should the sequential dose trajectory of the recruited
#' patients be plotted,  along with the observed toxicities? Defaults to FALSE.
#' @param file File name where the dose-response plots are stored,  in a pdf
#' format. The program will ammend the current sample size to the end of the
#' file name.
#' @param print.it Whether to print the graph to screen. Defaults to
#'   \code{print.it=TRUE}.
#' @param ... Further arguments passed to or from other methods
#' @return If \code{trajectory=FALSE}), the function invisibly returns list
#'   containing 2 objects of class `ggplot'.
#' @author Michael Sweeting \email{mjs212@@medschl.cam.ac.uk} (University of
#' Cambridge,  UK)
#' @seealso \code{\link{bcrm}}
#' @references Sweeting M.,  Mander A.,  Sabin T. \pkg{bcrm}: Bayesian Continual
#' Reassessment Method Designs for Phase I Dose-Finding Trials. \emph{Journal
#' of Statistical Software} (2013) 54: 1--26.
#' \url{http://www.jstatsoft.org/article/view/v054i13}
#' @method ggplot bcrm
#' @export
ggplot.bcrm <- function(data, mapping=aes(), file=NULL, each=FALSE, title="Posterior p(DLT) quantiles: 2.5%,  25%,  50%,  75%,  97.5% \n Diamond shows next recommended dose",
                      trajectory=FALSE, cols.traj=c(No="blue", Yes="orange"), print.it=TRUE, ..., environment=NULL){
  if (is.null(data$data$dose)){
    data$data$dose <- data$sdose
  }

  dose.label <- if(is.null(data$dose)) "Standardised dose" else "Dose"

  f <- which.f(data$ff)

  logit2 <- data$method %in% c("exact", "exact.sim") & data$ff == "logit2"

  if (trajectory){
    ggplot.bcrm_trajectory(data, cols=cols.traj)
  } else {
    # Plots when there is no loss vector
    a <- if(is.null(data$loss)){
      if(!each){
        if(logit2){
          ggplot.bcrm_logit2_noloss(data, dose.label=dose.label)
        } else {
          ggplot.bcrm_nologit2_noloss(data, dose.label=dose.label)
        }
      } else { # each is TRUE
        if(logit2){
          ggplot.bcrm_logit2_noloss_each(data, dose.label=dose.label)
        } else {
          ggplot.bcrm_nologit2_noloss_each(data, dose.label=dose.label)
        }
      }
    } else {
      if(!each){
        ggplot.bcrm_loss(data, dose.label=dose.label)
      } else {
        ggplot.bcrm_loss_each(data, dose.label=dose.label)
      }
    }

    # Get the barchart, if possible
    b <- try(ggplot.bcrm_barchart(data, cols=cols.barchart, dose.label=dose.label), silent=TRUE)
    if (class(b)[1] == "try-error"){
      b <- NULL
    }

    if (print.it){
      if(!each){
        grid::grid.newpage()
        grid::pushViewport(grid::viewport(layout=grid::grid.layout(2, 1)))
        vplayout <- function(x, y)  grid::viewport(layout.pos.row=x, layout.pos.col=y)
        print(a, vp=vplayout(1, 1))
        if(!is.null(b)) print(b, vp=vplayout(2, 1))
      } else {
        print(a)
      }
    }
    invisible(list(a, b))
  }
}

################################################################################
# Barchart

ggplot.bcrm_barchart <- function(data, mapping=aes(), cols=cols.barchart, dose.label){
  df2 <- data.frame(dose = factor(c(rep(data$dose, data$tox), rep(data$dose, data$notox)), levels=data$dose),
                    Outcome = factor(c(rep("DLT", sum(data$tox)), rep("No DLT", sum(data$notox))), levels=c("DLT", "No DLT")))

  p <-
  ggplot(data=df2, aes(x=dose, fill=Outcome)) +
    geom_bar() +
    xlab(dose.label) +
    ylab("Number") +
    scale_fill_hue(limits=c("DLT", "No DLT"))

  p
}

################################################################################
# No loss vector supplied

ggplot.bcrm_logit2_noloss <- function(data=NULL, mapping=aes(), dose.label){
  d <- data$data
  ndose <- d$ndose[[length(d$ndose)]]$ndose

  p <-
  ggplot(df, aes(x=dose, y=est)) +
    geom_point() +
    geom_hline(aes(yintercept=target.tox), data=df, col=4, linetype=2) +
    xlab(dose.label) +
    ylab("Probability of DLT") +
    ylim(0, 1) +
    ggtitle("Posterior point estimates \n Diamond shows next recommended dose") +
    geom_point(aes(x=dose, y=est), data=df[ndose, ], size=4, col=4, shape=9)

  p
}

ggplot.bcrm_nologit2_noloss <- function(data=NULL, mapping=aes(), dose.label){
  quants <- data$ndose[[length(data$ndose)]]$quantiles
  if (nrow(quants) != 5){
    stop("The quantiles in data$ndose should contain 5 rows")
  }

  df <- data.frame(dose = data$dose,
                   target.tox = data$target.tox,
                   est = data$ndose[[length(data$ndose)]]$est,
                   mean = data$ndose[[length(data$ndose)]]$mean,
                   q2.5 = quants["2.5%", ],
                   q25 = quants["25%", ],
                   q50 = quants["50%", ],
                   q75 = quants["75%", ],
                   q97.5 = quants["97.5%", ])

  p <-
  ggplot(df, aes(x=dose, ymin=q2.5, ymax=q97.5)) +
    geom_errorbar(colour="red") +
    geom_pointrange(aes(x=dose, y=q50, ymin=q25, ymax=q75), data=df, fill="red") +
    geom_hline(aes(yintercept=target.tox), data=df, col=4, linetype=2) +
    xlab(dose.label) +
    ylab("Probability of DLT") +
    ylim(0, 1) +
    ggtitle("Posterior p(DLT) quantiles: 2.5%,  25%,  50%,  75%,  97.5% \n Diamond shows next recommended dose") +
    geom_point(aes(x=dose, y=q50), data=df[data$ndose[[length(data$ndose)]][[1]], ], size=4, col=4, shape=9)

  p
}

################################################################################
# No loss vector supplied, each cohort

get.bcrm.data.noeach <- function(data){
  data.frame(dose = rep(data$dose, length(data$ndose)),
             target.tox = data$target.tox,
             cohort = rep(0:(length(data$ndose)-1), each=length(data$dose)),
             est = c(sapply(1:length(data$ndose), function(i){ data$ndose[[i]]$est })),
             ndose = rep(sapply(1:length(data$ndose), function(i){ dose[data$ndose[[i]]$ndose] }),
                         each=length(dose)))
}

ggplot.bcrm_logit2_noloss_each <- function(data, mapping=aes(), dose.label){
  df <- get.bcrm.data.noeach(data)

  p <-
  ggplot(df, aes(x=dose, y=est)) +
    geom_point()+
    geom_hline(aes(yintercept=target.tox), col=4, linetype=2)  +
    xlab(dose.label) +
    ylab("Probability of DLT") +
    ylim(0, 1) +
    ggtitle("Posterior point estimates \n Diamond shows next recommended dose") +
    geom_point(aes(x=ndose, y=q50), data=df[df$dose==df$ndose, ], size=4, col=4, shape=9) +
    facet_wrap(~ cohort)

  p
}

get.bcrm.data.each <- function(data){
  data.frame(dose = rep(data$dose, length(data$ndose)),
             target.tox = data$target.tox,
             cohort = rep(0:(length(data$ndose)-1), each=length(data$dose)),
             est = c(sapply(1:length(data$ndose), function(i){ data$ndose[[i]]$est })),
             mean = c(sapply(1:length(data$ndose), function(i){ data$ndose[[i]]$mean })),
             q2.5 = c(sapply(1:length(data$ndose), function(i){ data$ndose[[i]]$quantiles["2.5%", ] })),
             q25 = c(sapply(1:length(data$ndose), function(i){ data$ndose[[i]]$quantiles["25%", ] })),
             q50 = c(sapply(1:length(data$ndose), function(i){ data$ndose[[i]]$quantiles["50%", ] })),
             q75 = c(sapply(1:length(data$ndose), function(i){ data$ndose[[i]]$quantiles["75%", ] })),
             q97.5 = c(sapply(1:length(data$ndose), function(i){ data$ndose[[i]]$quantiles["97.5%", ] })),
             ndose = rep(sapply(1:length(data$ndose), function(i){ data$dose[data$ndose[[i]]$ndose] }), each=length(data$dose)))
}

ggplot.bcrm_nologit2_noloss_each <- function(data=NULL, mapping=aes(), dose.label){
  df <- get.bcrm.data.each(data)

  p <-
  ggplot(df, aes(x=dose, ymin=q2.5, ymax=q97.5)) +
    geom_errorbar(colour = "red") +
    geom_pointrange(aes(x = dose, y=q50, ymin=q25, ymax=q75), fill="red") +
    geom_hline(aes(yintercept = target.tox), data=df, col=4, linetype=2) +
    xlab(dose.label) +
    ylab("Probability of DLT") +
    ylim(0, 1) +
    ggtitle("Posterior p(DLT) quantiles: 2.5%,  25%,  50%,  75%,  97.5% \n Diamond shows next recommended dose") +
    geom_point(aes(x=ndose, y=q50), data=df[df$dose==df$ndose, ], size=4, col=4, shape=9) +
    facet_wrap(~ cohort)

  p
}

################################################################################
# Loss vector supplied

get.bcrm.interval.data <- function(data){
  data.frame(cohort=rep(0:(length(data$ndose) - 1), each=length(data$loss)),
             xmin=min(data$dose),
             xmax=max(data$dose),
             ymin=c(0, tox.cutpoints),
             ymax=c(data$tox.cutpoints, 1),
             Loss=data$loss)
}

ggplot.bcrm_loss <- function(data=NULL, mapping=aes(), dose.label){
  df <- get.bcrm.data.noeach(data)
  df.interval <- get.bcrm.interval.data(data)

  p <-
  ggplot(df, aes(x=dose, ymin=q2.5, ymax=q97.5), colour="red") +
    geom_errorbar() +
    geom_pointrange(aes(x=dose, y=q50, ymin=q25, ymax=q75), data=df, fill="red") +
    geom_point(aes(x=dose, y=q50), data=df[data$ndose[[length(data$ndose)]][[1]], ], size=4, col=4, shape=9) +
    geom_rect(aes(xmin=xmin, ymin=ymin, xmax=xmax, ymax=ymax, fill=Loss), data=df.intervals, alpha=0.3) +
    scale_fill_gradient(breaks=sort(unique(df.intervals$Loss)), high="red", low="#99ccff", guide="legend") +
    xlab(dose.label) +
    ylab("Probability of DLT") +
    ylim(0, 1) +
    ggtitle("Posterior p(DLT) quantiles: 2.5%,  25%,  50%,  75%,  97.5% \n Diamond shows next recommended dose")

  p
}

ggplot.bcrm_loss_each <- function(data=NULL, mapping=aes(), dose.label){
  df <- get.bcrm.data.each(data)
  df.interval <- get.bcrm.interval.data(data)

  p <-
  ggplot(df, aes(x=dose, ymin=q2.5, ymax=q97.5), colour="red") +
    geom_errorbar() +
    geom_pointrange(aes(x=dose, y=q50, ymin=q25, ymax=q75), data=df, fill="red") +
    geom_point(aes(x=ndose, y=q50), data=df[df$dose==df$ndose, ], size=4, col=4, shape=9) +
    geom_rect(aes(xmin=xmin, ymin=ymin, xmax=xmax, ymax=ymax, fill=Loss), data=df.intervals, alpha=0.3) +
    xlab(dose.label) +
    ylab("Probability of DLT") +
    ylim(0, 1) +
    ggtitle(title) +
    facet_wrap(~ cohort)

  p
}

################################################################################
# Trajectory

ggplot.bcrm_trajectory <- function(data=NULL, mapping=aes(), cols=c(No="blue", Yes="orange")){
  data <- data$data
  data$Toxicity <- ifelse(data$tox == 1, "Yes", "No")

  p <-
  ggplot(data, aes(x=patient, y=dose, shape=Toxicity, colour=Toxicity)) +
    geom_point(size=3) +
    scale_colour_manual(values=cols) +
    xlab("Patient") +
    ylab("Dose Level")

  p
}
