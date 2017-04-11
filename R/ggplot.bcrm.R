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
#' @param ... Further arguments passed to or from other methods
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
                      trajectory=FALSE, cols.traj=c(No="blue", Yes="orange"), ..., environment=NULL){
  dose <- if(is.null(data$dose)) data$sdose else data$dose
  dose.label <- if(is.null(data$dose)) "Standardised dose" else "Dose"
  f <- which.f(data$ff)

  if (trajectory){
    ggplot.bcrm.trajectory(data$data)
  }

  d <- get.bcrm.plot.data(data, each)
  df <- d[[1]]
  df2 <- d[[2]]

  a <- if(is.null(data$loss)){
    if(!each){
      if(data$method %in% c("exact", "exact.sim") & data$ff == "logit2"){
        ggplot(df, aes(x=dose, y=est)) +
          geom_point() +
          geom_hline(aes(yintercept=target.tox), data=df, col=4, linetype=2) +
          xlab(dose.label)+ylab("Probability of DLT") +
          ylim(0, 1) +
          ggtitle("Posterior point estimates \n Diamond shows next recommended dose") +
          geom_point(aes(x=dose, y=est), data=df[data$ndose[[length(data$ndose)]][[1]], ], size=4, col=4, shape=9)
      } else {
        ggplot(df, aes(x=dose, ymin=q2.5, ymax=q97.5)) +
          geom_errorbar(colour="red") +
          geom_pointrange(aes(x=dose, y=q50, ymin=q25, ymax=q75), data=df, fill="red") +
          geom_hline(aes(yintercept=target.tox), data=df, col=4, linetype=2) +
          xlab(dose.label) +
          ylab("Probability of DLT") +
          ylim(0, 1) +
          ggtitle("Posterior p(DLT) quantiles: 2.5%,  25%,  50%,  75%,  97.5% \n Diamond shows next recommended dose") +
          geom_point(aes(x=dose, y=q50), data=df[data$ndose[[length(data$ndose)]][[1]], ], size=4, col=4, shape=9)
      }
    } else { # each is TRUE
      if(data$method %in% c("exact", "exact.sim") & data$ff=="logit2"){
        ggplot()+geom_point(aes(x=dose, y=est), data=df)+
          geom_hline(aes(yintercept=target.tox), data=df, col=4, linetype=2)  +
          xlab(dose.label) +
          ylab("Probability of DLT") +
          ylim(0, 1) +
          ggtitle("Posterior point estimates \n Diamond shows next recommended dose") +
          geom_point(aes(x=ndose, y=q50), data=df[df$dose==df$ndose, ], size=4, col=4, shape=9) +
          facet_wrap(~ cohort)
      } else {
        ggplot() +
          geom_errorbar(aes(x=dose, ymin=q2.5, ymax=q97.5), colour="red", data=df) +
          geom_pointrange(aes(x=dose, y=q50, ymin=q25, ymax=q75), data=df, fill="red") +
          geom_hline(aes(yintercept=target.tox), data=df, col=4, linetype=2) +
          xlab(dose.label)+ylab("Probability of DLT") +
          ylim(0, 1) +
          ggtitle("Posterior p(DLT) quantiles: 2.5%,  25%,  50%,  75%,  97.5% \n Diamond shows next recommended dose") +
          geom_point(aes(x=ndose, y=q50), data=df[df$dose==df$ndose, ], size=4, col=4, shape=9) +
          facet_wrap(~ cohort)
      }
    }
  } else {
    df.intervals <- data.frame(cohort=rep(0:(length(data$ndose)-1), each=length(data$loss)),
                               xmin=min(dose),
                               xmax=max(dose),
                               ymin=c(0, tox.cutpoints),
                               ymax=c(data$tox.cutpoints, 1),
                               Loss=data$loss)
    if(!each){
      ggplot() +
        geom_errorbar(aes(x=dose, ymin=q2.5, ymax=q97.5), colour="red", data=df) +
        geom_pointrange(aes(x=dose, y=q50, ymin=q25, ymax=q75), data=df, fill="red") +
        geom_point(aes(x=dose, y=q50), data=df[data$ndose[[length(data$ndose)]][[1]], ], size=4, col=4, shape=9) +
        geom_rect(aes(xmin=xmin, ymin=ymin, xmax=xmax, ymax=ymax, fill=Loss), data=df.intervals, alpha=0.3) +
        scale_fill_gradient(breaks=sort(unique(df.intervals$Loss)), high="red", low="#99ccff", guide="legend") +
        xlab(dose.label) +
        ylab("Probability of DLT")+ylim(0, 1) +
        ggtitle("Posterior p(DLT) quantiles: 2.5%,  25%,  50%,  75%,  97.5% \n Diamond shows next recommended dose")
    } else {
      ggplot() +
        geom_errorbar(aes(x=dose, ymin=q2.5, ymax=q97.5), colour="red", data=df) +
        geom_pointrange(aes(x=dose, y=q50, ymin=q25, ymax=q75), data=df, fill="red") +
        geom_point(aes(x=ndose, y=q50), data=df[df$dose==df$ndose, ], size=4, col=4, shape=9) +
        geom_rect(aes(xmin=xmin, ymin=ymin, xmax=xmax, ymax=ymax, fill=Loss), data=df.intervals, alpha=0.3) +
        xlab(dose.label) +
        ylab("Probability of DLT")+ylim(0, 1) +
        ggtitle(title) +
        facet_wrap(~ cohort)
    }
  }
  b <- if(nrow(df2)!=0) {
    ggplot()+geom_bar(aes(x=dose, fill=Outcome), data=df2) +
      xlab(dose.label) +
      ylab("Number")+scale_fill_hue(limits=c("DLT", "No DLT"))
  } else { NULL  }

  if(!each){
    grid.newpage()
    pushViewport(viewport(layout=grid.layout(2, 1)))
    vplayout <- function(x, y)  viewport(layout.pos.row=x, layout.pos.col=y)
    print(a, vp=vplayout(1, 1))
    if(!is.null(b)) print(b, vp=vplayout(2, 1))
  } else {
    print(a)
  }
}

ggplot.bcrm.trajectory <- function(data=NULL, mapping=aes(), cols=c(No="blue", Yes="orange"), ..., environment=NULL){
  data$Toxicity <- ifelse(df$tox == 1, "Yes", "No")

  ggplot(data, aes(x=patient, y=dose, shape=Toxicity, colour=Toxicity)) +
    geom_point(size=3) +
    scale_colour_manual(values=cols) +
    xlab("Patient") +
    ylab("Dose Level")
}


get.bcrm.plot.data <- function(data, each){
  if(!each){
    if(data$method %in% c("exact", "exact.sim") & data$ff=="logit2"){
      df <- data.frame(dose = dose, target.tox = data$target.tox, est = data$ndose[[length(data$ndose)]]$est)
    } else {
      df <- data.frame(dose = dose,
                       target.tox = data$target.tox,
                       est = data$ndose[[length(data$ndose)]]$est,
                       mean = data$ndose[[length(data$ndose)]]$mean,
                       q2.5 = data$ndose[[length(data$ndose)]]$quantiles["2.5%", ],
                       q25 = data$ndose[[length(data$ndose)]]$quantiles["25%", ],
                       q50 = data$ndose[[length(data$ndose)]]$quantiles["50%", ],
                       q75 = data$ndose[[length(data$ndose)]]$quantiles["75%", ],
                       q97.5 = data$ndose[[length(data$ndose)]]$quantiles["97.5%", ])
    }
    df2 <- data.frame(dose = factor(c(rep(dose, data$tox), rep(dose, data$notox)), levels=dose),
                      Outcome = factor(c(rep("DLT", sum(data$tox)), rep("No DLT", sum(data$notox))), levels=c("DLT", "No DLT")))
  } else {
    if(data$method %in% c("exact", "exact.sim") & data$ff=="logit2"){
      df <- data.frame(dose = rep(dose, length(data$ndose)),
                       target.tox = data$target.tox,
                       cohort = rep(0:(length(data$ndose)-1), each=length(dose)),
                       est = c(sapply(1:length(data$ndose), function(i){data$ndose[[i]]$est})),
                       ndose = rep(sapply(1:length(data$ndose), function(i){dose[data$ndose[[i]]$ndose]}), each=length(dose)))
    } else {
      df <- data.frame(dose = rep(dose, length(data$ndose)),
                       target.tox = data$target.tox,
                       cohort = rep(0:(length(data$ndose)-1), each=length(dose)),
                       est = c(sapply(1:length(data$ndose), function(i){data$ndose[[i]]$est})),
                       mean = c(sapply(1:length(data$ndose), function(i){data$ndose[[i]]$mean})),
                       q2.5 = c(sapply(1:length(data$ndose), function(i){data$ndose[[i]]$quantiles["2.5%", ]})),
                       q25 = c(sapply(1:length(data$ndose), function(i){data$ndose[[i]]$quantiles["25%", ]})),
                       q50 = c(sapply(1:length(data$ndose), function(i){data$ndose[[i]]$quantiles["50%", ]})),
                       q75 = c(sapply(1:length(data$ndose), function(i){data$ndose[[i]]$quantiles["75%", ]})),
                       q97.5 = c(sapply(1:length(data$ndose), function(i){data$ndose[[i]]$quantiles["97.5%", ]})),
                       ndose = rep(sapply(1:length(data$ndose), function(i){dose[data$ndose[[i]]$ndose]}), each=length(dose)))
    }
    df2 <- data.frame()

    list(df, df2)
  }
}
