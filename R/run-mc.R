
#' Run a Monte Carlo simulation of a scenario
#' 
#' Run a Monte Carlo simulation for a scenario given fixed distributions for one
#' or more of the simulation parameters.  At each iteration all of the
#' parameters for which distributions were supplied will have values selected
#' from their distributions, while other parameters will be held at the base
#' values supplied.
#' 
#' The distparams input should be a named list of functions that return samples
#' from the indicated variables.  For example, if we wanted to have \code{hospFrac}
#' distributed as Beta(4, 98) and \code{T0} distributed as Gamma(20, 4), we would
#' write:
#' \preformatted{
#' distparms <- list(
#'    T0 = function(N) {rgamma(N, 20, 4)},
#'    hospFrac = function(N) {rbeta(N, 4, 98)})
#' }
#' Helper functions such as \code{\link{genbeta}} provide an easy way to generate
#' functions for common distributions.  Using the helper functions you would write
#' these distributions as:
#' \preformatted{
#' distparams <- list(
#'     T0 = gengamma(20, 4),
#'     hospFrac = genbeta(4, 98))
#' }
#' 
#' @param N number of MC iterations to run
#' @param baseparms Base parameters.  Any parameters without distributions will
#' be set to the values specified here (or to their default values, if unspecified.)
#' @param distparms Parameter distributions.  See details for how to specify a
#' distribution.
#' @param tmax End time for simulations
#' @param rawsims Flag indicating whether to return the raw simulation output, or
#' a summarized version.
#' @return A list of data frames, one for each of the forecast variables returned by 
#' \code{\link{run_scenario}.  The form of the outputs depends on the setting
#' of the \code{rawsims} input.  If \code{rawsims} is \code{FALSE}, then we get back
#' a data frame with ensemble mean, median, standard deviation, 5th, and 95th percentiles,
#' reported by time.
#' If \code{rawsims} is \code{TRUE}, then each output is a data frame with the output
#' variable reported by time and MC iteration.}
#' @importFrom foreach %do% %dopar%
#' @export
run_mc <- function(N, distparms, tmax, baseparms=list(), scenarioName='HospCensus', 
                   rawsims=FALSE)
{
  stopifnot(N>3)
  stopifnot(length(distparms) > 0)
  
  ## dummy variables to silence warnings
  scenario <- iter <- value <- NULL
  
  ## get the parameter values for the iterations
  parmvals <- lapply(distparms, function(f) {f(N)})
  runparms <- baseparms
  
  raw_rslts <-
    foreach::foreach(iter=seq(1,N), .combine=dplyr::bind_rows) %dopar% {
      parmdraws <- lapply(parmvals, function(l) {l[[iter]]})
      for(parm in names(parmdraws)) {
        runparms[[parm]] <- parmdraws[[parm]]
      }
      
      run_scenario(tmax, runparms, iter)
    }
  ## Scenario column currently contains the iteration number; rename and put
  ## the user's scenario name in the appropriate column.
  raw_rslts <- dplyr::rename(raw_rslts, mciter = scenario)
  raw_rslts$scenario <- scenarioName
  
  ## Set up some mappings between names of variables and columns in the results
  vars <- c('newCases', 'symptoInfection', 'acuteHosp', 'icuHosp', 'imvHosp')
  timevars <- list(newCases = 'time', symptoInfection = 'daysToSympto',
                  acuteHosp = 'daysToAcuteHosp', icuHosp = 'daysToicuHosp',
                  imvHosp = 'daysToIMVHosp')
  
  rslts <- 
    if(rawsims) {
      lapply(vars, function(v) {
        varout <- raw_rslts[c(timevars[[v]], v, 'mciter', 'scenario')]
        names(varout) <- c('time', 'value', 'mciter', 'scenario')
        varout[['varName']] <- v
        varout
      })
    } else {
      lapply(vars, function(v) {
        varraw <- raw_rslts[c(timevars[[v]], v, 'mciter', 'scenario')]
        names(varraw) <- c('time', 'value', 'mciter', 'scenario')
        dplyr::group_by(varraw, time, scenario) %>%
          dplyr::summarise(mean = mean(value),
                           median = median(value),
                           sdev = sd(value),
                           q05 = quantile(value, 0.05),
                           q95 = quantile(value, 0.95)) %>%
          dplyr::ungroup()
      })
    }
  names(rslts) <- vars
  rslts
}

#' Helper functions for creating parameter distributions
#' 
#' These functions create parameter distributions for use in \code{\link{run_mc}}.
#' 
#' @param alpha,beta Beta distribution parameters
#' @param shape,rate Gamma distribution parameters
#' @param lambda Poisson distribution parameter
#' @param mean,sd Normal distribution parameters
#' @name pdist
NULL

#' @describeIn pdist Create a beta distribution with specified parameters.
#' @export
genbeta <- function(alpha, beta) {function(N) {rbeta(N, alpha, beta)}}

#' @describeIn pdist Create a gamma distribution with specified shape and rate.
#' @export
gengamma <- function(shape, rate) {function(N) {rgamma(N, shape, rate)}}

#' @describeIn pdist Create a Poisson distribution with specified mean.
#' @export
genpois <- function(lambda) {function(N) {rpois(N, lambda)}}

#' @describeIn pdist Create a normal distribution with specified mean and standard deviation.
#' @export
gennorm <- function(mean, sd) {function(N) {rnorm(N, mean, sd)}}
