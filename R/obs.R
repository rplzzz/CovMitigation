#' Prepare the observed data for use in the likelihood function
#' 
#' @return List of observations used in the likelihood function.  This includes
#' \describe{
#'   \item{obsdata}{Confirmed new cases by date and county.  Also includes 
#'   an estimate of the number of people tested in each county each day.}
#'   \item{fips_codes}{Lookup table of FIPS codes for the Locality names output by
#'   the model.}
#'   \item{default_parm_vals}{Default values for the parameters used to compare
#'   model output to observed data (the parameters used inside the model have their
#'   own defaults)}
#' }
#' 
#' @param maxdate Maximum date to include in the comparison data.
#' @export
#' @keywords internal
get_obsdata <- function(maxdate = as.Date('2021-12-31'))
{
  ntesteff <- nposeff <- NULL
  obsdata <- vdhcovid::va_weekly_ntest_county
  obsdata$time <- as.numeric(obsdata$date - as.Date('2020-01-01'))
  obsdata <- dplyr::rename(obsdata, ntest=ntesteff, npos=nposeff)
  obsdata <- obsdata[obsdata$date <= maxdate, ]
  
  ## Ordinarily I would arrange the obs in the same order that the output comes
  ## from the model, but in this case we can't be guaranteed that we will have
  ## obs for every county at every model time, so we're just going to have to 
  ## do a join after each model run.
  
  ## Default values for likelihood parameters.  Ideally we should collect all of
  ## the parameter defaults together somewhere.
  default_parm_vals <- c(
    b0 = 10.0,
    b1 = 0.0
  )
  
  hosp <- hospitalizations

  list(obsdata=obsdata, hosp=hosp, default_parm_vals=default_parm_vals)
}


#' Create simulated observations from a scenario run
#' 
#' Given a matrix of simulation outputs and a function of time for the number
#' of tests, create simulated observations.  
#' 
#' @param rundata Matrix of run data returned from an integration of 
#' \code{\link{seir_equations}}
#' @param bparms Vector containing the b0 and b1 parameters.
#' @param population Total population for the simulated entity
#' @param ntestfn Function of time that returns the number of daily tests performed.
#' @export
simobs <- function(rundata, bparms, population=100000, ntestfn=linear_ntest)
{
  timecol <- which(colnames(rundata) == 'time')
  mdata <- rundata[ , -c(timecol)]
  totpop <- round(apply(mdata, 1, sum))
  totinfct <- round(mdata[,'I'] + mdata[,'Is'])
  fi <- totinfct / totpop
  
  time <- rundata[, timecol]
  ntest <- round(ntestfn(time))
  b <- bparms[['b0']] - bparms[['b1']]*log(ntest)
  biasedfi <- padjust(fi, b)
  
  N <- length(ntest)
  npos <- rbinom(N, ntest, biasedfi)
  
  d <- tibble::tibble(time=time, ntest=ntest, npos=npos)
  d$week <- ceiling(d$time / 7)
  d <- dplyr::group_by(d, week) %>%
    dplyr::summarise(time = max(time), ntest=sum(ntest), npos=sum(npos))
  
  d$locality <- 'simulated'  
  d$fips <- 99999
  d$population <- population
  
  d
}

## Linear ramp up of effective tests, starting 
linear_ntest <- function(t) {
  1 + t/10
}
