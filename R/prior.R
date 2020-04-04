
#' Generate a log prior function
#' 
#' This function returns a a function that given an input set of parameters evaluates
#' the log-prior probability density for the parameters.  Eventually this function
#' will be able to take hyperparameter arguments that influence the details of the
#' priors, but for now they're fixed.
#' 
#' The parameters currently recognized are described in \code{\link{gen_likelihood}}.
#' If there are any parameters not in that list present, an error will be raised.  
#' 
#' @return A function that takes a named vector of parameters and returns the log-prior
#' probability density
#' @export
gen_prior <- function()
{
  parm_names <- c('T0', 'D0', 'A0', 'Ts', 'day_zero', 'b')
  function(parms) {
    stopifnot(!is.null(names(parms)))
    if(any(!names(parms) %in% parm_names)) {
      badparms <- names(parms)[!names(parms) %in% parm_names]
      stop('unrecognized parameters: ', paste(badparms, collapse=', '))
    }
    ## For the time being we will want to run without Ts, so if it is missing,
    ## replace with the default value
    if(!'Ts' %in% names(parms)) {
      parms['Ts'] <- param_defaults[['Ts']]
    }
    
    ## All of our parameters have to be > 0, and b has to be < 1
    ## TODO make the sampler guarantee this.  How have we not included this in 
    ## the sampler already?
    if(any(parms <= 0)) {
      -Inf
    } else {
      logps <- c(
        dgamma(parms[c('T0', 'D0', 'A0', 'Ts')], 2, 0.5, log=TRUE),
        dnorm(parms['day_zero'], 60, 15, log=TRUE),
        dlnorm(parms['b'], 3, 2, log=TRUE)
      )
      sum(logps)
    }
  }
}