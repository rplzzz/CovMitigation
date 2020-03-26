
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
#' Optionally, instead of a number, beta may be a data frame with two columns,
#' 'time' and 'beta'.  The time column must start at zero and be strictly
#' increasing, while the beta column may hold any positive values. In this case,
#' beta is considered to be piecewise constant; each time t reaches the the next
#' value of 'time', the value of beta changes to the corresponding value from
#' the beta column.
#' 
#' @param t Simulation time value
#' @param variables Vector of current variable values (see details)
#' @param parameter List or vector of parameter values (see details)
#' @return A list, as described in \code{\link[deSolve]{ode}}.  In this case we provide only
#' the first element of the list, which is a vector of derivative values.
#' @export
sir_equations <- function(t, variables, parameters) {
  beta <- parameters[['beta']]
  gamma <- parameters[['gamma']]
  if(is.data.frame(beta)) {
    tmin <- min(beta$time)
    stopifnot(t >= tmin)
    irow <- rle(t >= beta$time)$lengths[1]    # Find the last beta$t that is <= t
    beta <- beta$beta[irow]
  }
  with(as.list(variables), {
    dS <- -beta * I * S
    dI <-  beta * I * S - gamma * I
    dR <-  gamma * I
    return(list(c(dS, dI, dR)))
  })
}
