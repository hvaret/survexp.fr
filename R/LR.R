#' Log-Rank test between an observed and an expected survival curve
#'
#' Log-Rank test between an observed and an expected survival curve
#'
#' @param futime follow-up time of the subjects in days
#' @param status 0 if censored or 1 if dead at \code{futime}
#' @param age age in days
#' @param sex \code{"male"} or \code{"female"}
#' @param entry_date entry date in the study
#' @param ratetable a table of event rates, such as \code{survexp.fr} or \code{survexp.us}
#' @return A list containing the observed number of deaths, the expected number of deaths, the Log-Rank statistic and its p-value
#' @author Hugo Varet
#' @details
#' The Log-Rank is calculated as:
#' \deqn{LR = (O-E)^2/E}
#' where \eqn{O} is the observed number of deaths and \eqn{E} is the expected number based on the patients' characteristics (sex, age and entry date in the study).
#' It follows a Khi-2 distribution with one degree of freedom, which allows to compute its p-value.
#' @references
#' R. Peto and J. Peto, Asymptotically Efficient Rank Invariant Test Procedures, Journal of the Royal Statistical Society, 1972
#' @examples
#' attach(data.example)
#' LR(futime, status, age, sex, entry_date)

LR=function(futime, status, age, sex, entry_date, ratetable=survexp.fr::survexp.fr){
  data <- na.omit(data.frame(futime, status, age, sex, entry_date))
  data$risk <- -log(survexp(futime~1,data=data,rmap=list(year=entry_date,age=age,sex=sex),cohort=F,ratetable=ratetable,conditional=TRUE))
  O <- sum(data$status)
  E <- sum(data$risk)
  LR <- ((O-E)^2)/E
  p.value <- 1-pchisq(LR,1)
  return(list(O=O, E=E, LR=LR, p.value=p.value))
}
