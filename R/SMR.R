#' Standardized Mortality Ratio (SMR)
#'
#' Computes the SMR, its confidence interval and its associated p-value
#'
#' @param futime follow-up time of the subjects in days
#' @param status 0 if censored or 1 if dead at \code{futime}
#' @param age age in days
#' @param sex \code{"male"} or \code{"female"}
#' @param entry_date entry date in the study
#' @param ratetable a table of event rates, such as \code{survexp.fr} or \code{survexp.us}
#' @param alpha determines the confidence level (1-\code{alpha}) of the confidence interval
#' @return A list containing the observed number of deaths, the expected number of deaths, the "classic" SMR
#' (with its confidence interval and its p-value) and the SMR calculated by a Poisson model (with its confidence interval and its p-value)
#' @author Jean-Philippe Jais and Hugo Varet
#' @details The SMR is estimated using two different methods.
#'
#' The classic method is:
#' \deqn{SMR = O/E}
#' where \eqn{O} is the observed number of deaths and \eqn{E} is the expected number based on the patients' characteristics (sex, age and entry date in the study).
#'
#' The SMR is also estimated performing a Poisson model where \eqn{O} is the dependant variable and \eqn{E} is an offset.
#'
#' @references
#' N. Breslow and N. Day, Statistical methods in cancer research, Volume II - The design and analysis of cohort studies, World Health Organization, 1987
#' @examples
#' attach(data.example)
#' SMR(futime, status, age, sex, entry_date)

SMR=function(futime, status, age, sex, entry_date,
             ratetable=survexp.fr::survexp.fr, alpha=0.05){
  data <- na.omit(data.frame(futime,status,age,sex,entry_date))
  data$risk <- -log(survexp(futime~1, data=data, 
                            rmap=list(year=entry_date, age=age, sex=sex),
                            cohort=FALSE, ratetable=ratetable,
                            conditional=TRUE))

  O <- sum(data$status)
  E <- sum(data$risk)

  # calcul classique
  SMR <- O/E
	SMR.lo <- O/E*(1-1/9/O-qnorm(1-alpha/2)/3/sqrt(O))^3
	SMR.up <- (O+1)/E*(1-1/9/(O+1)+qnorm(1-alpha/2)/3/sqrt(O+1))^3
  if (E>=10){
	  chisq <- ((abs(O-E)-0.5)^2)/E
  } else{
	  Mstar <- ifelse(O>=E, O, O+1)
	  chisq <- 9*Mstar*(1-(1/(9*Mstar))-((E/Mstar)^(1/3)))^2
	}
	p.value <- 1-pchisq(chisq,1)
  SMR.classic <- list(SMR=SMR, SMR.lo=SMR.lo, SMR.up=SMR.up, p.value=p.value)

  # Poisson model
  fit <- glm(O~1+offset(log(E)), family=poisson)
  coefs <- summary(fit)$coefficients
  SMR <- exp(coefs[1,1])
  SMR.lo <- exp(coefs[1,1] - qnorm(1-alpha/2)*coefs[1,2])
  SMR.up <- exp(coefs[1,1] + qnorm(1-alpha/2)*coefs[1,2])
  p.value <- coefs[1,4]
  SMR.poisson <- list(SMR=SMR,SMR.lo=SMR.lo,SMR.up=SMR.up,p.value=p.value)

  return(list(O=O, E=E, 
              SMR.classic=SMR.classic,
              SMR.poisson=SMR.poisson))
}
