
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
#' @param verbose If \code{TRUE}, print out additional diagnostic information.
#' @return A function that takes a named vector of parameters and returns the log-prior
#' probability density
#' @export
gen_prior <- function(hparms, verbose=FALSE)
{
  parm_names <- c('eta','xi', 'zeta', 'D0', 'A0', 'Ts', 'I0', 'mask_effect', 'b0', 'b1')
  function(parms) {
    stopifnot(!is.null(names(parms)))
    if(any(!names(parms) %in% parm_names)) {
      badparms <- names(parms)[!names(parms) %in% parm_names]
      stop('unrecognized parameters: ', paste(badparms, collapse=', '))
    }

    logps <- c(
      dgamma(parms[c('A0', 'Ts', 'D0')], 
             c(8,8, 14), 
             c(2,2, 2), log=TRUE),
      dnorm(parms['eta'], -0.7, 1, log=TRUE),
      dlnorm(parms['xi'], -0.5, 2, log=TRUE),
      dgamma(parms['zeta'], 1, 1, log=TRUE),         # mobility effect expected to be positive
      dlnorm(parms['I0'], 2, 1, log=TRUE),
      dgamma(-parms['mask_effect'], 1, 2, log=TRUE),  # mask effect expected to be negative
      dlnorm(parms['b0'], hparms[['bmulog']], hparms[['bsiglog']], log=TRUE),
      dexp(parms['b1'], 1, log=TRUE)
    )
    if(verbose) {
      pnames <- c('A0:\t', 'Ts:\t',  'D0:\t', 
                  'eta', 'xi', 'zeta',
                  'b:\t', 'I0:\t', 'mask_effect:\t')
      prvals <- signif(logps, 3)
      msg <- paste(pnames, prvals, collapse='\n')
      totmsg <- paste('total:\t', signif(sum(logps, na.rm=TRUE), 4))
      message('Prior contributions:\n', msg, '\n----------------\n', totmsg)
    }
    sum(logps, na.rm = TRUE)
  }
}

#' Evaluate the quantile for the multidimensional df defined by the priors
#' 
#' Given a vector of probabilities, find the corresponding quantile vector assuming
#' an independence copula and margins as described in \code{\link{gen_prior}}.
#' 
#' @section TO DO:
#' 
#' This function is currently inflexible with respect to the parameters included
#' in the calculation (because of the fixed list of parameter names).  Fix this
#' so that we can pass p as a named list.
#' 
#' @param p A vector of probabilities
#' @param hparms  Named vector of hyperparameters, described in \code{\link{gen_post}}
#' @export
qprior <- function(p, hparms=list()) {
  hparms <- fill_defaults(hparms, default_hparms)
  parm_names <- c('eta', 'xi', 'zeta', 'D0', 'A0', 'Ts', 'I0', 'mask_effect', 'b0','b1')
  stopifnot(length(p) == length(parm_names))
  
  ## Return a named vector of parameters
  c(
    eta = qnorm(p[1], -0.7, 1),
    xi = qlnorm(p[2], -0.5, 2),
    zeta = qgamma(p[3], 1, 1),
    D0 = qlnorm(p[4], 2, 0.5),
    A0 = qgamma(p[5], 8, 2),
    Ts = qgamma(p[6], 8, 2),
    I0 = qlnorm(p[7], 2, 1),
    mask_effect = -qgamma(p[8], 1, 2),
    b0 = qlnorm(p[9], hparms[['bmulog']], hparms[['bsiglog']]),
    b1 = qexp(p[10], 1)
  )
}
