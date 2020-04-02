
#' Prepare a likelihood function for use in Bayesian calibration
#' 
#' Build a function that given a set of model parameters runs the model and 
#' computes the log-likelihood for the model.
#' 
#' The observational data is the daily county-by-county reports of new confirmed
#' COVID-19 cases.  We are using the New York Times dataset
#' (\code{https://github.com/nytimes/covid-19-data}), which is based on reports
#' from state departments of health, including the VDH.  For statistics on COVID-19
#' testing, we use the \code{\link{uvads_covid19}} dataset, which was collected
#' directly from VDH by the UVA data science team.
#' 
#' The model parameters currently recognized are:
#' \describe{
#' \item{T0}{(real > 0) The initial doubling time for the number of infected people}
#' \item{D0}{(real > 0) The initial average infection duration}
#' \item{A0}{(real > 0) The average incubation time}
#' \item{Ts}{(real > 0) The average time from infection to symptom onset}
#' \item{day_zero}{(real > 0) The day on which the infection started in Fairfax County (not
#' necessarily when it was first observed), given in days since 1Jan2020.
#' This will be used to line up the observed dataset with model results.}
#' \item{b}{(real in (0,1)) The testing bias factor.  That is, 1/f, where f is the
#' factor by which testing results overestimate the infection fraction in the general
#' population.}
#' }
#' For the time being, we only accept a single \code{day_zero} parameter, and 
#' all of the counties are delayed relative to Fairfax by a number of days equal
#' to the difference between their first observed case and the first case observed
#' in Fairfax.
#' 
#' To compute the likelihood, our assumption is that the output of the model 
#' represents an average infection rate.  The average observation rate is then
#' \deqn{\bar{N}_o = \bar{N}_I f_t b,} where \eqn{f_t} is the fraction of the total
#' population tested.  We then assume that the actual observations are distributed
#' as \deqn{N_o \sim Pois(\bar{N}_o).}
#' 
#' 
#' @return A function that takes a named list of model parameters and returns
#' a log-likelihood value.  See details for the parameters recognized.
#' @importFrom dplyr %>%
#' @export
gen_likelihood <- function()
{
  ## silence warnings
  ntest <- ftest <- state <- date <- FIPS <- fips <- 
    cases <- day <- Locality <- locality <- newCases <- newcases <- NULL
  
  ## First get the datasets.  Filter and reformat as necessary.
  teststats <- 
    dplyr::filter(uvads_covid19, !is.na(ntest)) %>%
    dplyr::select(date, ftest)
  
  obsdata <- 
    dplyr::filter(NYTimesCOVID19::cov19county, 
                  state=='Virginia',
                  date >= min(teststats[['date']])) %>%
    dplyr::mutate(day = as.numeric(date - as.Date('2020-01-01')),
                  fips = as.integer(fips)) %>%
    dplyr::group_by(fips) %>%
    dplyr::mutate(newcases = c(0, diff(cases))) %>%    # add daily new cases
    dplyr::ungroup() %>%
    ## sometimes you get a decrease due to transfers or data corrections.  Treat
    ## these as zero
    dplyr::mutate(newcases = pmax(0, newcases)) %>%
    dplyr::select(date, fips, newcases, day)
    
  ## Our locality names are not the same as NYT, but we have a table of FIPS codes
  fips_codes <- dplyr::select(va_county_first_case, FIPS, Locality) %>%
    dplyr::rename(fips=FIPS)
  
  ## Ordinarily I would arrange the obs in the same order that the output comes
  ## from the model, but in this case we can't be guaranteed that we will have
  ## obs for every county at every model time, so we're just going to have to 
  ## do a join after each model run.
  
  function(parms) {
    ## The day-zero parameter is being handled in a simplified way compared to
    ## how it is handled in the model, so pull it out of the parameter list.
    modparms <- parms[names(parms) != 'day_zero']
    day0 <- parms[['day_zero']]
    
    ## model time corresponding to the observations.  Note these might not be integers.
    obsdata[['time']] <- obsdata[['day']] - day0
    ## get output for every day up to the last in the dataset.
    tmax <- max(obsdata[['time']])
    tvals <- c(0, seq(to = tmax, length.out = floor(tmax)))
    
    modout <- run_scenario(tvals, modparms) %>%
      dplyr::select(time, locality, newCases) %>%
      dplyr::rename(Locality = locality, model.newcases = newCases) %>%
      dplyr::left_join(fips_codes, by='Locality')
    
    cmp <- dplyr::full_join(obsdata, modout, by=c('time', 'fips'))
    
    ## If any rows have model.cases missing, it means that our day-zero put the start
    ## of the outbreak after the first observed case in that county.  Such parameter
    ## values have likelihood == 0.
    ## Conversely, if cases is missing, it means that we had a nonzero projection
    ## from the model on a day before we had any observed cases, which is perfectly
    ## fine, so fill in those cases with zeros.
    
    if(any(is.na(cmp[['model.cases']]))) {
      -Inf
    }
    else {
      miss <- is.na(cmp$newcases)
      cmp$newcases[miss] <- 0
      logl <- dpois(cmp$newcases, cmp$model.newcases, log=TRUE)
      
      sum(logl)
    }
  }
}