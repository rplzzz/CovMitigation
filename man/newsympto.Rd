% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim.R
\name{newsympto}
\alias{newsympto}
\title{Calculate incidence of symptomatic cases}
\usage{
newsympto(I, Is, R, dropfirst = FALSE)
}
\arguments{
\item{I}{Vector of asymptomatic infections over time.}

\item{Is}{Vector of symptomatic infections over time.}

\item{R}{Vector of recovered/removed population over time.}

\item{dropfirst}{If \code{TRUE}, drop the first entry (see details)}
}
\description{
Given a time series of symptomatic and asymptomatic infections, and recovered
population, calculate the new symptomatic cases at each time period.
}
\details{
New symptomatic cases are tricky because we have cases entering and leaving
the pool.  Basically, we look at the change in Is and back out the portion of
the change in R that is attributable to recovery from symptomatic cases.
Because we assume the recovery rates for symptomatic and asymptomatic cases
are equal, this is proportional to the ratio of symptomatic to asymptomatic
infections.

Formally:\cr
ras = I / Is\cr
deltaIs = newSympto - deltaRs\cr
deltaR  = deltaRs + deltaRa\cr
deltaRa / deltaRs = I/Is = ras\cr
   --> deltaR = (1+ras) deltaRs  --> deltaRs = deltaR/(1+ras)\cr
newSympto = deltaIs + deltaRs\cr
where deltaX is the change in X from the previous time step.

Because the results depend on differences from the previous time step, we
can't calculate the new cases for the first time step.  If the first time
step is the beginning of the infection, then we're looking at an initial
steady state, and the initial difference values are therefore zero, so that's
what we use.  If the first data point is in the middle of a run, then this
assumption is no good, and we have to discard the first value.  The
\code{dropfirst} argument allows this to happen automatically.
}
