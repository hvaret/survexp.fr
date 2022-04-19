\name{survexp_plot}
\alias{survexp_plot}
\title{Observed Kaplan-Meier, expected and relative survival curves}
\description{Displays the observed Kaplan-Meier, expected and relative survival curves}

\usage{
survexp_plot(futime, status, age, sex, entry_date, ratetable = survexp.fr, 
             main = "Observed and expected survival", 
             xlab = "Time (years)", ylab = "Survival", 
             col.km = "black", lwd.km = 2, lty.km = 1, conf.int.km = TRUE, 
             col.exp = "blue", lwd.exp = 2, lty.exp = 1, 
             main.rel = "Relative survival", ylab.rel = "Relative survival", 
             col.rel = "black", lwd.rel = 2, lty.rel = 1, 
             times = seq(0,max(futime, na.rm = TRUE)/365.241,length=6)[-1], alpha = 0.05,
             xscale = 365.241, ...)
}

\arguments{
  \item{futime}{follow-up time of the subjects in days}
  \item{status}{0 if censored or 1 if dead at \code{futime}}
  \item{age}{age in days}
  \item{sex}{\code{"male"} or \code{"female"}}
  \item{entry_date}{entry date in the study}
  \item{ratetable}{a table of event rates, such as \code{survexp.fr} or \code{survexp.us}}
  \item{main}{main title of the Kaplan-Meier and expected survivals plot}
  \item{xlab}{x-label of the plot}
  \item{ylab}{y-label of the plot}
  \item{col.km}{color of the observed survival curve}
  \item{lwd.km}{line width of the observed survival curve}
  \item{lty.km}{line type of the observed survival curve}
  \item{conf.int.km}{\code{TRUE} to display the confidence interval of the observed survival}
  \item{col.exp}{color of the expected survival curve}
  \item{lwd.exp}{line width of the expected survival curve}
  \item{lty.exp}{line type of the expected survival curve}
  \item{main.rel}{main title of the relative survival plot} 
  \item{ylab.rel}{y-label of the relative survival plot}  
  \item{col.rel}{color of the relative survival curve}
  \item{lwd.rel}{line width of the relative survival curve}
  \item{lty.rel}{line type of the relative survival curve}
  \item{times}{times to draw the confidence intervals of the relative survival}
  \item{alpha}{determines the confidence level (1-\code{alpha}) of the confidence intervals for the relative survival}
  \item{xscale}{see the \code{xscale} argument in \code{\link{plot.survfit}}}
  \item{\dots}{other arguments to be passed in \code{\link{plot.survfit}}}
}

\details{This function displays the observed and expected survivals, and the relative survival which is defined as:
\deqn{r(t) = exp(-exp(\beta) \times t)}
where \eqn{exp(\beta)} is the excess risk by time unit estimated by an additive Poisson model.
}

\value{A matrix containing the values of relative survivals and their confidence intervals for each time of \code{times}}

\references{
M. Pohar and J. Stare, Making relative survival analysis relatively easy, Computers in Biology and Medicine, 2007

M. Pohar and J. Stare, Relative survival analysis in R, Computers Methods and Programs in Biomedicine, 2006
}

\author{Hugo Varet}

\examples{
attach(data.example)
survexp_plot(futime, status, age, sex, entry_date)
}
