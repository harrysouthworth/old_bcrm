#' CRM INTERACTIVE. Conduct a CRM trial interactively allowing user to specify outcomes after each cohort has been recruited
#' @param tox number of successes (toxicities)
#' @param notox number of failed patients (no-toxicities)
#' @param ncurrent Current no. of patients in trial
#' @param cohort Cohort size
#' @param ndose Proposed dose level for next cohort
#' @param dose Dose labels for each dose
#' @export
crm.interactive <- function(tox, notox, ncurrent, cohort, ndose, dose){
  k  <-  length(tox)
  repeat {
    if(is.null(dose)){
      cat("\n\n RECOMMENDED DOSE LEVEL FOR PATIENTS ", ncurrent-cohort+1, " to ", ncurrent,  "IS:",  ndose[[1]])
      ans  <-  get.dose.level(k, ncurrent, cohort)
    } else {
      cat("\n\n RECOMMENDED DOSE FOR PATIENTS ", ncurrent-cohort+1, " to ", ncurrent,  "IS:",  dose[ndose[[1]]])
      ans  <-  get.dose(dose, ncurrent, cohort)
    }
    if (ans==-2)
      ans  <-  ifelse(is.null(dose), ndose[[1]], dose[ndose[[1]]])
    if (ans==0) {
      cat("\n\n EXIT AND RETURN THE RESULTS SO FAR?")
      cat("\n DO YOU REALLY WANT TO EXIT ? (Y/N)  ")
      yn  <-  readline()
      if (yn=="y" || yn=="Y") break
      else                    next
    }
    # y.j is toxicity information from the treatment of patient j in the current
    # cohort of patients at dose level ndose
    y <- vector()
    cat("\n")
    for(j in (ncurrent-cohort+1):ncurrent){
      cat("ENTER TOXICITY DATA FOR PATIENT", j, "(1=TOX,  0=NO TOX): ")
      y  <-  c(y, get.answer())
    }
    # give the user a last chance to modify the treatment and outcome
    cat("\n\n\t\t ENTERED VALUES:")
    if(is.null(dose)){
      cat("\n DOSE LEVEL ...",  ans)
    } else {
      cat("\n DOSE ...",  ans)
    }
    cat("\n TOXICITIES ....", y)
    cat("\n PRESS `RETURN' IF OK,  OR ANY OTHER KEY TO ENTER NEW VALUES ")
    key  <-  readline()
    if (nchar(key)==0) break
  }
  if (ans==0)  return(list(bk=TRUE))
  if(!is.null(dose)){
    ans <- which(dose==ans)
  }
  notox[ans] <- notox[ans]+(cohort-sum(y))
  tox[ans] <- tox[ans]+sum(y)
  return(list(tox=tox, notox=notox, y=y, bk=FALSE))
}
