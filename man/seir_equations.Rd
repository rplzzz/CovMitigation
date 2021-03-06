% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim.R
\name{seir_equations}
\alias{seir_equations}
\title{Differential equations for the SEIR model}
\usage{
seir_equations(t, variables, parameters)
}
\arguments{
\item{t}{Simulation time value}

\item{variables}{Vector of current variable values (see details)}

\item{parameters}{List or vector of parameter values (see details)}
}
\value{
A list, as described in \code{\link[deSolve]{ode}}.  In this case we provide only
the first element of the list, which is a vector of derivative values.
}
\description{
Computes the derivatives of the susceptible (S), exposed (E), infected (I),
symptomatic (Is), and recovered (R)
populations.  This version of the model omits the population birth and death rates.
}
\details{
The variables for this model are S, E, I, Is, and R, the susceptible, infected
asymptomatic, infected symptomatic, and recovered
populations.  They should be passed as a named vector, in that order.  (The order
is important because the ODE solver ignores names.)

The parameters for the model are:
\describe{
  \item{alpha}{Progression rate parameter (i.e., the rate at which people move from exposed to infected)}
  \item{beta}{Infection rate parameter:  relative rate at which infected people infect susceptible people.}
  \item{gamma}{Recovery rate parameter:  relative rate at which infected people recover (whether symptomatic or not)}
  \item{epsilon}{Symptom rate parameter:  relative rate at which asymptomatic people develop symptoms.}
}
These should be passed in as a named vector; the order doesn't matter.  Also,
note that previous versions absorbed the 1/N factor into beta, but this one
does not.

Optionally, instead of a number, alpha, beta and/or gamma may be data frames with two columns,
'time' and 'value'.  The time column must start at zero and be strictly
increasing, while the beta column may hold any positive values. In this case,
the parameter is considered to be piecewise constant; each time t reaches the the next
value of 'time', the value of the parameter changes to the corresponding value from
the value column.  (Time varying epsilon is not currently supported.)
}
