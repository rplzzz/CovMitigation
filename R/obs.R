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
#' @keywords internal
get_obsdata <- function()
{
  ntesteff <- nposeff <- NULL
  obsdata <- vdhcovid::va_weekly_ntest_county
  obsdata$time <- as.numeric(obsdata$date - as.Date('2020-01-01'))
  obsdata <- dplyr::rename(obsdata, ntest=ntesteff, npos=nposeff)
  
  ## Ordinarily I would arrange the obs in the same order that the output comes
  ## from the model, but in this case we can't be guaranteed that we will have
  ## obs for every county at every model time, so we're just going to have to 
  ## do a join after each model run.
  
  ## Default values for likelihood parameters.  Ideally we should collect all of
  ## the parameter defaults together somewhere.
  default_parm_vals <- c(
    b = 10.0
  )
  
  list(obsdata=obsdata, default_parm_vals=default_parm_vals)
}
