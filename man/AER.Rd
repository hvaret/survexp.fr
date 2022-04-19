\name{AER}
\alias{AER}
\title{Absolute Excess Risk (AER)}
\description{Computes the AER, its confidence interval and its associated p-value}
\usage{
AER(futime, status, age, sex, entry_date, PY.stand = 10000, 
    ratetable = survexp.fr, alpha = 0.05)
}

\arguments{
  \item{futime}{follow-up time of the subjects in days}
  \item{status}{0 if censored or 1 if dead at \code{futime}}
  \item{age}{age in days}
  \item{sex}{\code{"male"} or \code{"female"}}
  \item{entry_date}{entry date in the study}
  \item{PY.stand}{value to get the AER for \code{stand} person-years}
  \item{ratetable}{a table of event rates, such as \code{\link{survexp.fr}} or \code{\link{survexp.us}}}
  \item{alpha}{determines the confidence level (1-\code{alpha}) of the confidence interval}
}

\details{The Absolute Excess Risk (AER) is defined as:
\deqn{AER = O-E}
where \eqn{O} is the observed number of deaths and \eqn{E} is the expected number based on the patients'characteristics (sex, age and entry date in the study). This function uses an
additive Poisson model to compute the AER.
}

\value{A list containing the AER with the corresponding number of person-years (\code{PY.stand} argument), its confidence interval, its p-value, 
the observed number of deaths, the expected number of deaths and the observed number of person-years}

\references{
N. Breslow and N. Day, Statistical methods in cancer research, Volume II - The design and analysis of cohort studies, World Health Organization, 1987

P. Dickman, A. Sloggett, M. Hills and T. Hakulinen, Regression models for relative survival, Statistics in Medicine, 2004

C. Elie, Y. De Rycke, J.-P. Jais and P. Landais, Appraising relative and excess mortality in population-based studies of chronic diseases such as end-stage renal disease, Clinical Epidemiology, 2011
}

\author{Jean-Philippe Jais and Hugo Varet}

\examples{
attach(data.example)
AER(futime, status, age, sex, entry_date)
}
