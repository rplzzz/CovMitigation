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
  ## silence warnings
  ntest <- ntest <- state <- date <- FIPS <- fips <- 
    cases <- time <- Locality <- newcases <- NULL
  
  teststats <- dplyr::select(va_county_ntest, date, fips, ntest)
  
  obsdata <- 
    dplyr::filter(NYTimesCOVID19::cov19county, 
                  state=='Virginia',
                  !is.na(fips),
                  date >= min(teststats[['date']])) %>%
    dplyr::mutate(time = as.numeric(date - as.Date('2020-01-01')),
                  fips = as.integer(fips)) %>%
    dplyr::group_by(fips) %>%
    dplyr::mutate(newcases = c(0, diff(cases))) %>%    # add daily new cases
    dplyr::ungroup() %>%
    ## sometimes you get a decrease due to transfers or data corrections.  Treat
    ## these as zero
    dplyr::mutate(newcases = pmax(0, newcases)) %>%
    ## Add the estimates of number of tests
    dplyr::left_join(teststats, by=c('date', 'fips')) %>%
    dplyr::select(date, fips, newcases, ntest, time)
    
  ## Our locality names are not the same as NYT, but we have a table of FIPS codes
  fips_codes <- dplyr::select(va_county_first_case, FIPS, Locality) %>%
    dplyr::rename(fips=FIPS)
  
  ## Ordinarily I would arrange the obs in the same order that the output comes
  ## from the model, but in this case we can't be guaranteed that we will have
  ## obs for every county at every model time, so we're just going to have to 
  ## do a join after each model run.
  
  ## Default values for likelihood parameters.  Ideally we should collect all of
  ## the parameter defaults together somewhere.
  
  default_parm_vals <- c(
    day_zero = 60,
    b = 0.5
  )
  
  list(obsdata=obsdata, fips_codes=fips_codes, 
       default_parm_vals=default_parm_vals)
}
