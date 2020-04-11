
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
#' \item{Ts}{(real > 0) The average time from infection to symptom onset.  Note that
#' this parameter is probably not identifiable with our current data.}
#' \item{day_zero}{(real > 0) The day on which the infection started in Fairfax County (not
#' necessarily when it was first observed), given in days since 1Jan2020.
#' This will be used to line up the observed dataset with model results.}
#' \item{b}{(real > 0) The testing bias factor.  That is, the positive test rate
#' divided by the true infection rate.  b can be different from 1 because of false positives
#' or because testing is targeted to people suspected of having the disease.  (Generally,
#' we expect b>1, but we don't require this.)}
#' \item{I0}{(real > 0) The initial number of infected people, once the infection starts.
#' This is taken to be the same in all counties.}
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
#' @param fixed_parms A named vector of parameters to use as fixed values for parameters
#' not passed to the likelihood function.  Note these "fixed" parameters can still
#' be overridden by passing an explicit value.
#' @return A function that takes a named list of model parameters and returns
#' a log-likelihood value.  See details for the parameters recognized.
#' @importFrom dplyr %>%
#' @export
gen_likelihood <- function(fixed_parms=NULL)
{
  obs <- get_obsdata()
  obsdata <- obs$obsdata
  fips_codes <- obs$fips_codes

  if(!is.null(fixed_parms)) {
    for(p in names(fixed_parms)) {
      obs$default_parm_vals[[p]] <- fixed_parms[[p]]
    }
  }
  default_parm_vals <- obs$default_parm_vals
  
  function(parms) {
   ## complete the parameters from the defaults
    parms <- fill_defaults(parms, default_parm_vals)
    
    ## Separate out the parameters that get passed to the SEIR model.
    modparms <- as.list(parms[! names(parms) %in% c('day_zero', 'b')])
    if('day_zero' %in% names(parms)) {
      day0 <- parms[['day_zero']]
    }
    else {
      day0 <- default_parm_vals[['day_zero']]
    }
    
    ## get output for every day up to the last in the dataset.
    tmax <- max(obsdata[['time']])
    
    tvals <- c(day0, seq(ceiling(day0), tmax))
    
    ## Run the model and do a little post-processing
    modout <- run_scenario(tvals, modparms)
    mdata <- modout[c('time', 'locality','I','Is','population')]
    mdata$Itot <- mdata$I + mdata$Is
    mdata$fi <- mdata$Itot / mdata$population
    mdata <- dplyr::left_join(mdata, obs$fips_codes, by=c(locality='Locality'))
    
    cmp <- dplyr::full_join(obs$obsdata, mdata, by=c('time', 'fips'))
    cmp <- cmp[!(is.na(cmp$fi) | is.na(cmp$ntest)),]
    cmp <- cmp[(cmp$fi > 0) & (cmp$ntest > 0),]

    ## If any rows have Itot missing, it means that our day-zero put the start
    ## of the outbreak after the first observed case in that county.  Such parameter
    ## values have likelihood == 0.
    ## Conversely, if cases is missing, it means that we had a nonzero projection
    ## from the model on a day before we had any observed cases, which is perfectly
    ## fine, so fill in those cases with zeros.
    
    if(any(is.na(cmp$Itot))) {
      -Inf
    }
    else {
      miss <- is.na(cmp$newcases)
      cmp$newcases[miss] <- 0
      
      ## adjust model outputs for testing fraction and testing bias!
      b <- parms[['b']]

      ## The model forecast for the number of new cases is the current infection
      ## fraction (fi) times the number of tests performed.  However, we think that
      ## testing is biased toward people who are infected, so we multiply the odds
      ## ratio by a bias factor b
      biased_odds <- cmp$fi/(1-cmp$fi) * b
      cmp$biased_fi <- biased_odds / (1 + biased_odds)
      total_pop_va <- vaMyAgeBands$Total[1]
      popfac <- 1/total_pop_va
      cmp$popfrac <- cmp$population * popfac
      
      dhargs <- tibble::tibble(x=cmp$newcases,
                               m=cmp$biased_fi * cmp$population,
                               n=(1-cmp$biased_fi) * cmp$population,
                               k=(cmp$ntest * cmp$popfrac))
      
      ## There are a few counties and times where the number of observed cases is greater than
      ## the number of tests we're estimating.  In those cases, set the number of tests
      ## to the number of cases.
      dhargs$k <- pmax(dhargs$k, dhargs$x)
      
      ## Hypergeometric distribution.  Assume that tests performed in each county
      ## are proportional to the county's share of the total population.
      logl <- dhyper(dhargs$x,
                     ceiling(dhargs$m),
                     ceiling(dhargs$n),
                     ceiling(dhargs$k),
                     log=TRUE)
      
      if(any(is.na(logl))) {
        ## This seems to be happening occasionally.  Not sure why
        bad <- cmp[is.na(logl),]
        warning('Bad values in likelihood function.')
        for (i in seq(1,nrow(bad))) {
          warning(paste(colnames(bad), collapse='\t'), '\n', 
                  paste(bad[i,], collapse='\t'))
        }
      }
      sum(logl, na.rm=TRUE)
    }
  }
}

#' Prepare a posterior log-pdf function for use in Bayesian calibration
#' 
#' 
#' @param prior_weight Factor to multiply the prior by.  Default is to make it equal
#' to the number of counties for which we have confirmed COVID-19 cases.
#' @param fixed_parms A named vector of parameters to use as fixed values for parameters
#' not passed to the likelihood function.  Note these "fixed" parameters can still
#' be overridden by passing an explicit value.
#' @export
gen_post <- function(prior_weight=NULL, fixed_parms=NULL)
{
  lprior <- gen_prior()
  llik <- gen_likelihood(fixed_parms)
  
  if(is.null(prior_weight)) {
    prior_weight <- sum(!is.na(va_county_first_case$firstDay))
  }
  
  function(parms) {
    logp <- prior_weight * lprior(parms)
    if(is.finite(logp)) {
      logp <- logp + llik(parms)
    }
    logp
  }
}


#' Prepare the observed data for use in the likelihood function
#' 
#' @return List of observations used in the likelihood function.  This includes
#' \describe{
#'   \item{obsdata}{Confirmed new cases by date and county.  Also includes 
#'   statewide fraction of people tested.}
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
  
  ## First get the datasets.  Filter and reformat as necessary.
  teststats <- 
    dplyr::filter(uvads_covid19, !is.na(ntest)) %>%
    dplyr::select(date, ntest)
  
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
    dplyr::left_join(teststats, by='date') %>%
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