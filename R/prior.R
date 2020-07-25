
## distribution parameters for priors
hprior <- list(
    D0gshp = 8, D0grate = 2,
    A0gshp = 8, A0grate = 2,
    Tsgshp = 8, Tsgrate = 2,
    I0lmu = 2, I0lsig = 1,
    MEgshp = 1, MEgrate = 2,
    b0lmu = 4, b0lsig = 0.5,
    b1exprate = 1,

    betalmu = -0.7, betalsig = 1
)

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
  common_parm_names <- c('D0', 'A0', 'Ts', 'I0', 'mask_effect', 'b0', 'b1')
  function(parms) {
    stopifnot(!is.null(names(parms)))

    lbeta <- extract_beta_parms(names(parms), 'logical')
    beta_parms <- parms[lbeta]
    common_parms <- parms[!lbeta]

    if(any(!names(common_parms) %in% common_parm_names)) {
      badparms <- names(common_parms)[!names(common_parms) %in% common_parm_names]
      stop('unrecognized parameters: ', paste(badparms, collapse=', '))
    }

    logps <- c(
      dlnorm(beta_parms, hprior$betalmu, hprior$betalsig, log=TRUE),
      dgamma(common_parms['A0'], hprior$A0gshp, hprior$A0grate, log=TRUE),
      dgamma(common_parms['Ts'], hprior$Tsgshp, hprior$Tsgrate, log=TRUE),
      dgamma(common_parms['D0'], hprior$D0gshp, hprior$D0grate, log=TRUE),
      dlnorm(common_parms['I0'], hprior$I0lmu, hprior$I0lsig, log=TRUE),
      dgamma(-common_parms['mask_effect'], hprior$MEgshp, hprior$MEgrate, log=TRUE),  # mask effect expected to be negative
      dlnorm(common_parms['b0'], hprior$b0lmu, hprior$b0lsig, log=TRUE),
      dexp(common_parms['b1'], hprior$b1exprate, log=TRUE)
    )
    if(verbose) {
        pnames <- c(paste0(names(beta_parms),'\t'),
                    'A0:\t', 'Ts:\t',  'D0:\t',
                    'I0:\t', 'mask_effect:\t',
                    'b0:\t', 'b1:\t')
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
  common_parm_names <- c('D0', 'A0', 'Ts', 'I0', 'mask_effect', 'b0','b1')

  lbeta <- extract_beta_parms(names(p), 'logical')
  beta_p <- p[lbeta]
  common_p <- p[!lbeta]

  stopifnot(length(common_p) == length(common_parm_names))

  beta_vals <- qlnorm(beta_p, hprior$betalmu, hprior$betalsig)
  names(beta_vals) <- names(p)[lbeta]

  ## Return a named vector of parameters
  c(
    beta_vals,
    D0 = qgamma(common_p[['D0']], hprior$D0gshp, hprior$D0grate),
    A0 = qgamma(common_p[['A0']], hprior$A0gshp, hprior$A0grate),
    Ts = qgamma(common_p[['Ts']], hprior$Tsgshp, hprior$Tsgrate),
    I0 = qlnorm(common_p[['I0']], hprior$I0lmu, hprior$I0lsig),
    mask_effect = -qgamma(common_p[['mask_effect']], hprior$MEgshp, hprior$MEgrate),
    b0 = qlnorm(common_p[['b0']], hprior$b0lmu, hprior$b0lsig),
    b1 = qexp(common_p[['b1']], hprior$b1exprate)
  )
}
