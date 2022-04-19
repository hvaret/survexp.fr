\name{SMR}
\alias{SMR}
\title{Standardized Mortality Ratio (SMR)}
\description{Computes the SMR, its confidence interval and its associated p-value}
\usage{
SMR(futime, status, age, sex, entry_date, ratetable = survexp.fr, alpha = 0.05)
}

\arguments{
  \item{futime}{follow-up time of the subjects in days}
  \item{status}{0 if censored or 1 if dead at \code{futime}}
  \item{age}{age in days}
  \item{sex}{\code{"male"} or \code{"female"}}
  \item{entry_date}{entry date in the study}
  \item{ratetable}{a table of event rates, such as \code{survexp.fr} or \code{survexp.us}}
  \item{alpha}{determines the confidence level (1-\code{alpha}) of the confidence interval}
}

\details{The SMR is estimated using two different methods.

The classic method is:
\deqn{SMR = O/E}
where \eqn{O} is the observed number of deaths and \eqn{E} is the expected number based on the patients' characteristics (sex, age and entry date in the study).

The SMR is also estimated performing a Poisson model where \eqn{O} is the dependant variable and \eqn{E} is an offset.
}

\value{A list containing the observed number of deaths, the expected number of deaths, the "classic" SMR (with its confidence interval and its p-value) and the SMR calculated by a Poisson model (with its confidence interval and its p-value)}

\references{N. Breslow and N. Day, Statistical methods in cancer research, Volume II - The design and analysis of cohort studies, World Health Organization, 1987}

\author{Jean-Philippe Jais and Hugo Varet}

\examples{
attach(data.example)
SMR(futime, status, age, sex, entry_date)
}
