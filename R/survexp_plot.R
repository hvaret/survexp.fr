#' Observed Kaplan-Meier, expected and relative survival curves
#'
#' Displays the observed Kaplan-Meier, expected and relative survival curves
#'
#' @param futime follow-up time of the subjects in days
#' @param status 0 if censored or 1 if dead at \code{futime}
#' @param age age in days
#' @param sex \code{"male"} or \code{"female"}
#' @param entry_date entry date in the study
#' @param ratetable a table of event rates, such as \code{survexp.fr} or \code{survexp.us}
#' @param main main title of the Kaplan-Meier and expected survivals plot
#' @param xlab x-label of the plot
#' @param ylab y-label of the plot
#' @param col.km color of the observed survival curve
#' @param lwd.km line width of the observed survival curve
#' @param lty.km line type of the observed survival curve
#' @param conf.int.km \code{TRUE} to display the confidence interval of the observed survival
#' @param col.exp color of the expected survival curve
#' @param lwd.exp line width of the expected survival curve
#' @param lty.exp line type of the expected survival curve
#' @param main.rel main title of the relative survival plot
#' @param ylab.rel y-label of the relative survival plot
#' @param col.rel color of the relative survival curve
#' @param lwd.rel line width of the relative survival curve
#' @param lty.rel line type of the relative survival curve
#' @param times times to draw the confidence intervals of the relative survival
#' @param alpha determines the confidence level (1-\code{alpha}) of the confidence intervals for the relative survival
#' @param xscale see the \code{xscale} argument in \code{\link{plot.survfit}}
#' @param \dots other arguments to be passed in \code{\link{plot.survfit}}
#' @return A matrix containing the values of relative survivals and their confidence intervals for each time of \code{times}
#' @author Hugo Varet
#' @details This function displays the observed and expected survivals, and the relative survival which is defined as:
#' \deqn{r(t) = exp(-exp(\beta) \times t)}
#' where \eqn{exp(\beta)} is the excess risk by time unit estimated by an additive Poisson model.
#' @references
#' M. Pohar and J. Stare, Making relative survival analysis relatively easy, Computers in Biology and Medicine, 2007
#'
#' M. Pohar and J. Stare, Relative survival analysis in R, Computers Methods and Programs in Biomedicine, 2006
#' @examples
#' attach(data.example)
#' survexp_plot(futime, status, age, sex, entry_date)

survexp_plot=function(futime, status, age, sex, entry_date, 
                      ratetable=survexp.fr::survexp.fr,
                      main="Observed and expected survival", xlab="Time (years)", ylab="Survival",
                      col.km="black", lwd.km=2, lty.km=1, conf.int.km=TRUE,
                      col.exp="blue", lwd.exp=2, lty.exp=1,
                      main.rel="Relative survival", ylab.rel="Relative survival",
                      col.rel="black", lwd.rel=2, lty.rel=1, 
                      times=seq(0, max(futime,na.rm=TRUE)/365.241, length=6)[-1],
                      alpha=0.05, xscale=365.241, ...){
  data <- na.omit(data.frame(futime, status, age, sex, entry_date))
  par(mfrow=c(1,2))
  # observed survival
  km <- survfit(Surv(futime,status)~1,data=data)
  plot(km, main=main, xlab=xlab, ylab=ylab, lwd=lwd.km, col=col.km,
       conf.int=conf.int.km, xscale=xscale, ...)
  # expected survival
  expected <- survexp(futime~1, rmap=list(sex=sex, year=entry_date, age=age),
                      conditional=TRUE, data=data, ratetable=ratetable)
  lines(expected, lwd=lwd.exp, col=col.exp, lty=lty.exp, xscale=xscale)
  # legend
  legend("bottomleft",
         col=c(col.km, col.exp),
         lty=c(lty.km, lty.exp),
         lwd=c(lwd.km, lwd.exp),
         legend=c("Observed", "Expected"))
  # relative
  fit <- AER(data$futime, data$status, data$age, data$sex, data$entry_date,
             PY.stand=1, ratetable=ratetable, alpha=alpha)
  coef <- log(fit$AER)
  coef.lo <- log(fit$AER.lo)
  coef.up <- log(fit$AER.up)
  surv.rel <- exp(-exp(coef)*times)
  surv.rel.lo <- exp(-exp(coef.lo)*times)
  surv.rel.up <- exp(-exp(coef.up)*times)
  fun_curve <- function(x){exp(-exp(coef)*x)}
  curve(fun_curve, from=0, to=1.1*max(times), main=main.rel, xlab=xlab, ylab=ylab.rel,
        col=col.rel, lwd=lwd.rel, lty=lty.rel, ylim=c(0,1))
  arrows(times, surv.rel.lo, times, surv.rel.up, col=col.rel,
         lwd=lwd.rel, angle=90, code=3, length=0.07)
  points(times, surv.rel, pch=19, col=col.rel)
  return(cbind(times, surv.rel, surv.rel.lo, surv.rel.up))
}