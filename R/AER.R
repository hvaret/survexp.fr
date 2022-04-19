#' Absolute Excess Risk (AER)
#'
#' Computes the AER, its confidence interval and its associated p-value
#'
#' @param futime follow-up time of the subjects in days
#' @param status 0 if censored or 1 if dead at \code{futime}
#' @param age age in days
#' @param sex \code{"male"} or \code{"female"}
#' @param entry_date entry date in the study
#' @param PY.stand value to get the AER for \code{stand} person-years
#' @param ratetable a table of event rates, such as \code{survexp.fr} or \code{survexp.us}
#' @param alpha determines the confidence level (1-\code{alpha}) of the confidence interval
#' @return A list containing the AER with the corresponding number of person-years (\code{PY.stand} argument), its confidence interval, its p-value,
#' the observed number of deaths, the expected number of deaths and the observed number of person-years
#' @author Jean-Philippe Jais and Hugo Varet
#' @details The Absolute Excess Risk (AER) is defined as:
#' \deqn{AER = O-E}
#' where \eqn{O} is the observed number of deaths and \eqn{E} is the expected number based on the patients'characteristics (sex, age and entry date in the study).
#' This function uses an additive Poisson model to compute the AER.
#' @references
#' N. Breslow and N. Day, Statistical methods in cancer research, Volume II - The design and analysis of cohort studies, World Health Organization, 1987
#'
#' P. Dickman, A. Sloggett, M. Hills and T. Hakulinen, Regression models for relative survival, Statistics in Medicine, 2004
#'
#' C. Elie, Y. De Rycke, J.-P. Jais and P. Landais, Appraising relative and excess mortality in population-based studies of chronic diseases such as end-stage renal disease, Clinical Epidemiology, 2011
#' @examples
#' attach(data.example)
#' AER(futime, status, age, sex, entry_date)

AER=function(futime, status, age, sex, entry_date, PY.stand=10000,
             ratetable=survexp.fr::survexp.fr, alpha=0.05){
  data <- na.omit(data.frame(futime, status, age, sex, entry_date, 
                             person_year=futime/365.241))
  data$risk <- -log(survexp(futime~1, data=data, 
                            rmap=list(year=entry_date,age=age,sex=sex),
                            cohort=FALSE, ratetable=ratetable,
                            conditional=TRUE))
  mypoisson <- poisson()
  mypoisson$link <- "glm with Poisson error and non standard link function"
  mypoisson$linkfun <- function(mu) log(mu-E)
  mypoisson$linkinv <- function(eta) exp(eta)+E
  mypoisson$initialize <- expression({
    if (any(y < 0)) stop("Negative values not allowed for the Poisson family")
    n <- rep.int(1, nobs)
    mustart <- pmax(y, linkinv(-1000)) + 0.1                                    # linkinv(-1000) returns expected
  })

  E <- sum(data$risk)
  O <- sum(data$status)
  PY.observed <- sum(data$person_year)

  fit <- glm(O~1+offset(log(PY.observed)), family=mypoisson)
  coefs <- summary(fit)$coefficients
  AER <- exp(coefs[1,1]+log(PY.observed))*PY.stand/PY.observed                     # le +log(PY) et /PY s'annulent
  AER.lo <- exp(coefs[1,1]-qnorm(1-alpha/2)*coefs[1,2]+log(PY.observed))*PY.stand/PY.observed
  AER.up <- exp(coefs[1,1]+qnorm(1-alpha/2)*coefs[1,2]+log(PY.observed))*PY.stand/PY.observed
  return(list(AER=AER, PY.stand=PY.stand, 
              AER.lo=AER.lo, AER.up=AER.up, 
              p.value=coefs[1,4], O=O, E=E, 
              PY.observed=PY.observed))
}
