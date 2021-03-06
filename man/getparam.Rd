% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim.R
\name{getparam}
\alias{getparam}
\title{Get a time variable parameter from a table of step-changes}
\usage{
getparam(t, param)
}
\arguments{
\item{t}{Simulation time value}

\item{param}{A data frame giving the step changes in the parameter over time, or a single
constant parameter value}
}
\value{
The current value for the parameter, if time variable, or the constant value, if not.
}
\description{
Check to see if a parameter was given as a table of changes over time.  If so,
get the value of the parameter for the current time.  Otherwise, return the parameter's
constant value.
}
\keyword{internal}
