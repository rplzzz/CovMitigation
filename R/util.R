#' Fill in default parameter values
#' 
#' Works on vectors and lists.  Returns a structure of the same type as its first input.
#' 
#' @param parms Named vector of parameters
#' @param defaults Named vector of default parameters
#' @export
fill_defaults <- function(parms, defaults)
{
  for(n in names(defaults)) {
    if(! n %in% names(parms)) {
      parms[[n]] <- defaults[[n]]
    }
  }
  parms
}

#' Caclulate model coefficients from parameters
#' 
#' We parameterize the model SEIR in terms of initial doubling time, incubation period, 
#' and recovery time.  The model equations are given in terms of rates,
#' \eqn{\alpha}, \eqn{\beta}, and \eqn{\gamma}.  Also of interest is the basic
#' reproduction parameter \eqn{R_0}.
#' 
#' Recovery rate and progression rate are just the reciprocals of their respective
#' periods:
#' \deqn{\alpha = 1/A_0}
#' \deqn{\gamma = 1/D_0}
#' 
#' Beta is calculated from \eqn{\alpha}, \eqn{\gamma}, and \eqn{t_d} as
#' \deqn{\beta = \gamma + \ln 2 \frac{(\alpha+\gamma)}{\alpha t_d}}
#' 
#' \eqn{R_0} is calculated from \eqn{\beta} and \eqn{\gamma}.  Note that in the 
#' absence of mortality during the exposed period, \eqn{R_0} is independent of 
#' \eqn{\alpha}.
#' \deqn{R_0 = \frac{\beta}{\gamma}}
#' 
#' @param parms Named vector of model parameters
#' @param alpha SEIR alpha parameter
#' @param beta SEIR beta parameter
#' @param gamma SEIR gamma parameter
#' @name coefs
NULL


#' @describeIn coefs Calculate progression rate \eqn{\alpha} from model parameters
#' 
#' @export
calcalpha <- function(parms)
{
  1/parms[['A0']]
}

#' @describeIn coefs Calculate transmissibility \eqn{\beta} from model parameters.
#' 
#' This is the approximate version, using the Taylor series expansion.
#' 
#' @export
calcbeta_approx <- function(parms)
{
  alpha <- calcalpha(parms)
  gamma <- calcgamma(parms)
  td <- parms[['T0']]
  gamma + log(2) *(alpha+gamma) / (alpha*td)
}

#' @describeIn coefs Calculate recovery rate \eqn{\gamma} from model parameters.
#' 
#' @export
calcgamma <- function(parms)
{
  1/parms[['D0']]
}

#' @describeIn coefs Calculate beta, given model input parameters
#' 
#' This is the exact version, solving the nonlinear equation for beta.
#' 
#' @export
calcbeta <- function(parms)
{
  beta_guess <- calcbeta_approx(parms)
  alpha <- calcalpha(parms)
  gamma <- calcgamma(parms)
  targfun <- function(beta) {
    parms[['T0']] - calct0(alpha, beta, gamma)
  }
  solv <- nleqslv::nleqslv(beta_guess, targfun)
  if(solv$termcd != 1) {
    warning('calcbeta: failed to achieve convergence: msg= ', solv$message)
  }
  solv$x
}

#' @describeIn coefs Calculate basic reproduction number \eqn{R_0} from model parameters.
#' 
#' @export
calcr0 <- function(parms)
{
  calcbeta(parms)/calcgamma(parms)
}

#' @describeIn coefs Calculate the growth rate (lplus) given alpha, beta, gamma coefficients
#' 
#' @export
calclplus <- function(alpha, beta, gamma)
{
  apg <- alpha+gamma
  0.5*(-apg + sqrt(apg^2 + 4*alpha*(beta-gamma)))
}

#' @describeIn coefs Calculate the doubling time (t0), given alpha, beta, gamma coefficients
#' 
#' @export
calct0 <- function(alpha, beta, gamma)
{
  log(2) / calclplus(alpha, beta, gamma)
}