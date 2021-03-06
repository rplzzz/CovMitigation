% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bayes-filter.R
\name{run_parmset}
\alias{run_parmset}
\title{Run the compartment model for a single set of parameters}
\usage{
run_parmset(parms, istate, timevals, localityscen = NULL)
}
\arguments{
\item{parms}{Named vector of model parameters: beta, A0, D0, Ts, and mask_effect}

\item{istate}{Named vector of initial values for state variable: S, E, I, Is, and R}

\item{timevals}{Vector of times to output results.  The values of the state
variables will be initialized at the first time in the list.}

\item{localityscen}{\code{LocalScen} object returned from
\code{\link{localize_scenario}}.  If \code{NULL}, then no changes to
parameters will be applied in the future scenario.}
}
\value{
Matrix with a column for time and a column for each state variable.
The first row of the matrix (which would contain the initial time) is dropped
so that rbinding the output of several consecutive calls to this function will
not produce any duplicate rows.
}
\description{
By contrast to \code{\link{run_parms}}, this function runs just the ODE
integration for the SEIR equations.
}
