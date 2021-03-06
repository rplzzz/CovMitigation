% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/likelihood.R
\name{calc_nhosp}
\alias{calc_nhosp}
\title{Calculate the expected number of hospitalizations for a model}
\usage{
calc_nhosp(ph, modout, mfadjust = TRUE)
}
\arguments{
\item{ph}{Per-person hospitalization probability}

\item{modout}{Model output from \code{\link{run_scenario}}}

\item{mfadjust}{If \code{TRUE}, adjust for market fraction}
}
\value{
Data frame of time and expected hospitalizations
}
\description{
Calculate the expected number of hospitalizations by day, based on an assumed
fraction of symptomatic patients that eventually go to the hospital.
}
\details{
To make this calculation, we have to convert this probability to a rate.  If
the recovery rate is \eqn{\gamma = 1/D0}, then the hospitalization rate \eqn{r_h} that gives
the desired probability \eqn{P_h} is
\deqn{r_h = \gamma \frac{P_h}{1-P_h}}
}
