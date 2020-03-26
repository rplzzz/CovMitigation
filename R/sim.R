
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
#' Optionally, instead of a number, beta and gamma may be data frames with two columns,
#' 'time' and 'value'.  The time column must start at zero and be strictly
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
sir_equations <- function(t, variables, parameters)
{
  beta <- getparam(t, parameters[['beta']])
  gamma <- getparam(t, parameters[['gamma']])
  with(as.list(variables), {
    dS <- -beta * I * S
    dI <-  beta * I * S - gamma * I
    dR <-  gamma * I
    return(list(c(dS, dI, dR)))
  })
}

#' Get a time variable parameter from a table of step-changes
#' 
#' Check to see if a parameter was given as a table of changes over time.  If so,
#' get the value of the parameter for the current time.  Otherwise, return the parameter's
#' constant value.
#' 
#' @param t Simulation time value
#' @param param A data frame giving the step changes in the parameter over time, or a single
#' constant parameter value
#' @return The current value for the parameter, if time variable, or the constant value, if not.
#' @keywords internal
getparam <- function(t, param)
{
  if(is.data.frame(param)) {
    tmin <- min(param$time)
    stopifnot(t >= tmin)
    irow <- rle(t >= param$time)$lengths[1]    # Find the last param$t that is <= t
    param$value[irow]
  }
  else {
    param
  }
}
