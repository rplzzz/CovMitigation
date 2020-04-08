
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
  default_parm_vals <- obs$default_parm_vals
  if(!is.null(fixed_parms)) {
    for(p in names(fixed_parms)) {
      default_parm_vals[[p]] <- fixed_parms[[p]]
    }
  }
  
  function(parms) {
   ## complete the parameters from the defaults
    pnames <- names(default_parm_vals)
    use_defaults <- pnames[!pnames %in% names(parms)]
    for(n in use_defaults) {
      parms[[n]] <- default_parm_vals[[n]]      
    }
    
    ## Separate out the parameters that get passed to the SEIR model.
    modparms <- as.list(parms[! names(parms) %in% c('day_zero', 'b')])
    if('day_zero' %in% names(parms)) {
      day0 <- parms[['day_zero']]
    }
    else {
      day0 <- default_parm_vals[['day_zero']]
    }
    
    ## model time corresponding to the observations.  Note these might not be integers.
    obsdata[['time']] <- obsdata[['day']] - day0
    ## get output for every day up to the last in the dataset.
    tmax <- max(obsdata[['time']])
    tvals <- c(0, seq(to = tmax, length.out = floor(tmax)))
    
    modout <- run_scenario(tvals, modparms)
    cmp <- align_modout(modout, list(obsdata=obsdata, fips_codes=fips_codes, 
                                     default_parm_vals=default_parm_vals))

    ## If any rows have model.cases missing, it means that our day-zero put the start
    ## of the outbreak after the first observed case in that county.  Such parameter
    ## values have likelihood == 0.
    ## Conversely, if cases is missing, it means that we had a nonzero projection
    ## from the model on a day before we had any observed cases, which is perfectly
    ## fine, so fill in those cases with zeros.
    
    if(any(is.na(cmp$model.newcases))) {
      -Inf
    }
    else {
      miss <- is.na(cmp$newcases)
      cmp$newcases[miss] <- 0
      
      ## adjust model outputs for testing fraction and testing bias!
      if('b' %in% names(parms)) {
        b <- parms[['b']]
      }
      else {
        b <- default_parm_vals[['b']]
      }
      cmp$model.newcases <- cmp$model.newcases * cmp$ftest * b
      ## Occasionally the model will produce very small predictions that can
      ## cause problems in dpois.  Set a floor under the model output that is
      ## small enough that essentially any observations will be a no-go
      cmp$model.newcases <- pmax(cmp$model.newcases, 1e-8)
      
      logl <- dpois(cmp$newcases, cmp$model.newcases, log=TRUE)
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
#' @param fixed_parms A named vector of parameters to use as fixed values for parameters
#' not passed to the likelihood function.  Note these "fixed" parameters can still
#' be overridden by passing an explicit value.
#' @export
gen_post <- function(fixed_parms=NULL)
{
  lprior <- gen_prior()
  llik <- gen_likelihood(fixed_parms)
  function(parms) {
    logp <- lprior(parms)
    if(is.finite(logp)) {
      logp <- logp + llik(parms)
    }
    logp
  }
}

#' Align model output to observed data for comparison
#' 
#' @param modout Raw model output
#' @param obs Observed data, as returned by \code{\link{get_obsdata}}.  The time
#' column (= day - day0) must already have been added.
#' @keywords internal
align_modout <- function(modout, obs)
{
  mdata <- modout[c('time','locality','newCases')]
  names(mdata) <- c('time','Locality','model.newcases')
  mdata <- dplyr::left_join(mdata, obs$fips_codes, by='Locality')
  mdata$time <- round(mdata$time, 2)  ## fix roundoff error
    
  obs$obsdata$time <- round(obs$obsdata$time, 2)
  cmp <- dplyr::full_join(obs$obsdata, mdata, by=c('time', 'fips')) 
  
  cmp[!is.na(cmp$ftest),]
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
  ntest <- ftest <- state <- date <- FIPS <- fips <- 
    cases <- day <- Locality <- newcases <- NULL
  
  ## First get the datasets.  Filter and reformat as necessary.
  teststats <- 
    dplyr::filter(uvads_covid19, !is.na(ntest)) %>%
    dplyr::select(date, ftest)
  
  obsdata <- 
    dplyr::filter(NYTimesCOVID19::cov19county, 
                  state=='Virginia',
                  !is.na(fips),
                  date >= min(teststats[['date']])) %>%
    dplyr::mutate(day = as.numeric(date - as.Date('2020-01-01')),
                  fips = as.integer(fips)) %>%
    dplyr::group_by(fips) %>%
    dplyr::mutate(newcases = c(0, diff(cases))) %>%    # add daily new cases
    dplyr::ungroup() %>%
    ## sometimes you get a decrease due to transfers or data corrections.  Treat
    ## these as zero
    dplyr::mutate(newcases = pmax(0, newcases)) %>%
    dplyr::left_join(teststats, by='date') %>%
    dplyr::select(date, fips, newcases, ftest, day)
    
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