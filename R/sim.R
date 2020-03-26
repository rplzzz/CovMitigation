
#' Differential equations for the SIR model
#' 
#' Computes the derivatives of the susceptible (S), infected (I), and recovered (R) 
#' populations.
#' 
#' The variables for this model are S, I, and R, the susceptible, infected, and recovered 
#' populations.  They should be passed as a named vector, in that order.
#' 
#' The parameters for the model are beta, the infection parameter, and gamma, the recovery
#' parameter.  These should be passed in as a named vector; the order doesn't matter.
#' 
#' @param time Simulation time value
#' @param variables Vector of current variable values (see details)
#' @param parameter Vector of parameter values (see details)
#' @return A list, as described in \code{\link[deSolve]{ode}}.  In this case we provide only
#' the first element of the list, which is a vector of derivative values.
#' @export
sir_equations <- function(time, variables, parameters) {
  with(as.list(c(variables, parameters)), {
    dS <- -beta * I * S
    dI <-  beta * I * S - gamma * I
    dR <-  gamma * I
    return(list(c(dS, dI, dR)))
  })
}
