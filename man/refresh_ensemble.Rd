% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/refresh-ensemble.R
\name{refresh_ensemble}
\alias{refresh_ensemble}
\title{Refresh a model ensemble by replacing low performing members.}
\usage{
refresh_ensemble(
  modlist,
  mcmcdir,
  frac = 0.05,
  replacement_ensemble_override = NULL
)
}
\arguments{
\item{modlist}{List of filter model fits by county.}

\item{mcmcdir}{Name of the directory containing the MCMC inputs}

\item{frac}{Fraction of ensemble members to refresh (default = 0.05)}

\item{replacement_ensemble_override}{List of replacement ensembles. This is a
debugging parameter and will be removed.}
}
\description{
Identify the lowest performing members of an ensemble and replace them with 
random draws from the original MCMC.
}
