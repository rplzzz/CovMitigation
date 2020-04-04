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
#' @export
calcbeta <- function(parms)
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

#' @describeIn coefs Calculate basic reproduction number \eqn{R_0} from model parameters.
#' 
#' @export
calcr0 <- function(parms)
{
  calcbeta(parms)/calcgamma(parms)
}
