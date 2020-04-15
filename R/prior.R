
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
#' @param hparms Named vector of hyperparameters, described in \code{\link{gen_post}}
#' @return A function that takes a named vector of parameters and returns the log-prior
#' probability density
#' @export
gen_prior <- function(hparms)
{
  parm_names <- c('T0_hi', 'T0', 'D0', 'A0', 'Ts', 'day_zero', 'b', 'I0')
  function(parms) {
    stopifnot(!is.null(names(parms)))
    if(any(!names(parms) %in% parm_names)) {
      badparms <- names(parms)[!names(parms) %in% parm_names]
      stop('unrecognized parameters: ', paste(badparms, collapse=', '))
    }

    ## All of our parameters have to be > 0
    if(any(parms <= 0)) {
      -Inf
    } else {
      logps <- c(
        dgamma(parms[c('T0_hi', 'T0', 'A0', 'Ts')], 
               c(rep(hparms[['t0shp']],2), 8,8), 
               c(rep(hparms[['t0rate']],2), 2,2), log=TRUE),
        dlnorm(parms['D0'], 2, 0.5, log=TRUE),
        dnorm(parms['day_zero'], 30, 30, log=TRUE),
        dlnorm(parms['b'], hparms[['bmulog']], hparms[['bsiglog']], log=TRUE),
        dlnorm(parms['I0'], 2, 1, log=TRUE)
      )
      sum(logps, na.rm = TRUE)
    }
  }
}
