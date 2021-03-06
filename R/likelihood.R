
#' Prepare a likelihood function for use in Bayesian calibration
#'
#' Build a function that given a set of model parameters runs the model and
#' computes the log-likelihood for the model.
#'
#' The observational data is the daily county-by-county reports of new confirmed
#' COVID-19 cases.  We are using the New York Times dataset
#' (\code{https://github.com/nytimes/covid-19-data}), which is based on reports
#' from state departments of health, including the VDH.  For statistics on COVID-19
#' testing, we use the \code{\link[vacovdata]{vacovdata}} dataset, which was compiled
#' by the COVID Tracking project (\code{https://covidtracking.com/}).
#'
#' The model parameters currently recognized are:
#' \describe{
#' \item{eta}{(real) Log of the baseline transmissibility}
#' \item{xi}{(real) Coefficient of population density in log-transmissibility}
#' \item{D0}{(real > 0) The initial average infection duration}
#' \item{A0}{(real > 0) The average incubation time}
#' \item{Ts}{(real > 0) The average time from infection to symptom onset.  Note that
#' this parameter is probably not identifiable with our current data.}
#' \item{b}{(real > 0) The testing bias factor.  That is, the positive test rate
#' divided by the true infection rate.  b can be different from 1 because of false positives
#' or because testing is targeted to people suspected of having the disease.  (Generally,
#' we expect b>1, but we don't require this.)}
#' \item{I0}{(real > 0) The initial number of infected people, once the infection starts.
#' This is taken to be the same in all counties.}
#' }
#' For the time being, we start tracking the infection in Fairfax county on day 30, and
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
#' @section Likelihood function:
#'
#' By default, the function produced will return a single value, which is the sum of
#' \itemize{
#' \item{The log-likelihood function, summed over all data points}
#' \item{A correction for the hospitalization fraction}
#' \item{A correction for the symptomatic fraction}
#' }
#' Technically, the latter two are priors, since they don't depend on the observations,
#' but the model must be run in order to compute them, so they are computed here.
#'
#' If the \code{waicmode} flag is set, then the function will return detailed information
#' on the contribution of each data point to the log-likelihood.  The result will be a
#' data frame containing
#' \itemize{
#' \item{date}
#' \item{locality}
#' \item{expected counts}
#' \item{observed counts}
#' \item{Hypergeometric parameters x, m, n, and k}
#' \item{log likelihood}
#' }
#' Additionally, the antibody prevalence and symptomatic fraction corrections
#' will be attached to the data frame as attributes (\code{hfadjust} and
#' \code{fsadjust}).
#'
#' The default mode output can be used in optimizations or Markov chain Monte Carlo.
#' The waicmode output is useful for computing the WAIC, or for diagnosing which
#' counties and/or days are favoring one model over another.
#'
#' @param hparms Hyperparameters for the calculation, as described in \code{\link{gen_post}}
#' @param fixed_parms A named vector of parameters to use as fixed values for parameters
#' not passed to the likelihood function.  Note these "fixed" parameters can still
#' be overridden by passing an explicit value.
#' @param maxdate Last date to use in the calibration.  Default is 2020-06-30
#' @param verbose If \code{TRUE}, print additional diagnostic information.
#' @param waicmode If \code{TRUE}, generate a function that returns the likelihood
#' contributions from each data point.
#' @return A function that takes a named list of model parameters and returns
#' a log-likelihood value.  See details for the parameters recognized and the
#' output of the function returned.
#' @importFrom dplyr %>%
#' @export
gen_likelihood <- function(hparms, fixed_parms=NULL, maxdate=NULL, verbose=FALSE, waicmode=FALSE)
{
  ## silence warnings
  fips <- ntest <- biased_fi <- locality <- npos <- expected <-
    x <- m <- n <- k <- NULL

  if(is.null(maxdate)) {
    maxdate <- as.Date('2020-06-30')
  }
  obs <- get_obsdata(maxdate = maxdate)
  obsdata <- obs$obsdata

  if(!is.null(fixed_parms)) {
    for(p in names(fixed_parms)) {
      obs$default_parm_vals[[p]] <- fixed_parms[[p]]
    }
  }
  default_parm_vals <- obs$default_parm_vals

  ## Final dates of the weeks we use for aggregation
  wkdates <- sort(unique(obs$obsdata$date))
  ## Table of fips codes
  fips_table <- unique(dplyr::select(obs$obsdata, locality, fips))
  jan01 <- as.Date('2020-01-01')

  tantibody <- hparms[['tantibody']]
  fantibody <- hparms[['fantibody']]
  sigantibody <- hparms[['sigantibody']]

  function(parms) {
    ## complete the parameters from the defaults
    parms <- fill_defaults(parms, default_parm_vals)

    ## Separate out the parameters that get passed to the SEIR model.
    seirparms <- as.list(parms[! names(parms) %in% c('b0', 'b1')])

    ## get output for every day up to the last in the dataset.  Also, we need
    ## output at least up through the end of July to apply the antibody
    ## prevalence constraint
    tmax <- pmax(max(obsdata[['time']]), tantibody)
    t1 <- infection_t0             # constant defined in sim.R
    tvals <- seq(t1, tmax)

    ## Run the model and do a little post-processing
    modout <- run_scenario(tvals, seirparms)
    if(nrow(modout)==0) {
      ## All of the model outputs were after the the observation period.  Parameters
      ## are clearly bogus.
      return(-Inf)
    }
    fs <- fsympto(modout)
    mdata <- modout[c('time', 'locality','S', 'I','Is','population')]
    mdata$Itot <- mdata$I + mdata$Is
    mdata$fi <- mdata$Itot / mdata$population
    mdata <- dplyr::left_join(mdata, fips_table, by='locality')
    mdata$date <- mdata$time + jan01

    ## Compute the adjustment for antibody fraction
    antibody_data <- mdata[mdata$time == tantibody, c('locality', 'S', 'population')]
    fa <- 1 - antibody_data$S / antibody_data$population
    fa <- pmin(1, pmax(0, fa))          # fa sometimes gets out of the 0-1
                                        # interval by about 1e-8.
    correction <- losig(fa, fantibody, sigantibody)
    antibody_adjust <- sum(correction)

    if(verbose) {
      msg <- paste('\t',antibody_data$locality, ':\t',
                   signif(fa,3), '\t', signif(correction,4),
                   collapse='\n')
      message('Antibody prevalence adjustments:\n \tLocality \t\t fa \t correction\n',
              msg)
    }

    ## Aggregate the data by week.  For model outputs we want the average over
    ## the week.
    mdatawk <-
      dplyr::group_by(mdata, locality, fips) %>%
      dplyr::group_map(
        function(df, group) {
          wkagg(wkdates, df, c('I','Is','It'))
      },
      keep=TRUE) %>%
      dplyr::bind_rows()

    cmp <- dplyr::full_join(obs$obsdata, mdatawk, by=c('time','date', 'locality', 'fips'))
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
      miss <- is.na(cmp$npos)
      cmp$npos[miss] <- 0

      ## adjust model outputs for testing fraction and testing bias!
      ## Assume that the enrichment factor is never less than 1 (i.e., we never
      ## select *against* infection)
      b <- pmax(parms[['b0']] - parms[['b1']] * log(cmp$ntest), 1)

      ## The model forecast for the number of new cases is the current infection
      ## fraction (fi) times the number of tests performed.  However, we think that
      ## testing is biased toward people who are infected, so we multiply the odds
      ## ratio by a bias factor b
      cmp$biased_fi <- padjust(cmp$fi, b)

      dhargs <- tibble::tibble(x=cmp$npos,
                               m=cmp$biased_fi * cmp$population,
                               n=(1-cmp$biased_fi) * cmp$population,
                               k=cmp$ntest)

      ## Sometimes you get a number of observations greater than the total number
      ## of cases the model is projecting.  In theory, that model should be excluded
      ## with a logl of -Inf, but that can make it hard to find a good initial guess.
      ## Instead, adjust x so that that x <= m, but apply a finite penalty for each x > m.
      dhpenalty <- -1000 * sum(dhargs$x > ceiling(dhargs$m))
      dhargs$x <- pmin(dhargs$x, ceiling(dhargs$m))

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

      ## Compute a likelihood term for the statewide symptomatic fraction, which
      ## should be somewhere near 0.5.  Since this is just one value for the whole
      ## dataset, we weight it by the number of times in the data.  Note this term
      ## really only constrains Ts, and it is pretty much the only thing that constrains
      ## that parameter.
      fsadjust <-
        (max(obsdata$time) - min(obsdata$time)) *
        dbeta(fs, hparms[['fsalpha']], hparms[['fsbeta']], log=TRUE)

      if(verbose) {
        suml <- signif(sum(logl), 5)
        fsa <- signif(fsadjust, 5)
        msg <- paste(c('logl:\t', 'fsadjust:\t', 'antibody_adjust:\t', 'dhpenalty:\t'),
                     c(suml, fsa, antibody_adjust, dhpenalty), collapse='\n')
        message('Likelihood contributions:\n', msg)
      }

      if(waicmode) {
        rslt <- dplyr::bind_cols(cmp, dhargs) %>%
          dplyr::mutate(expected=biased_fi*ntest, logl=logl) %>%
          dplyr::select(date, fips, locality, npos, ntest, expected,
                        x, m, n, k, logl)
        attr(rslt, 'fsadjust') <- fsadjust
        attr(rslt, 'antibody_adjust') <- antibody_adjust
        rslt
      }
      else {
        sum(logl, na.rm=TRUE) + fsadjust + dhpenalty + antibody_adjust
      }
    }
  }
}

#' Prepare a posterior log-pdf function for use in Bayesian calibration
#'
#' Create a function that computes and sums the log-prior and the log-posterior.
#'
#' Recognized hyperparameters are
#' \describe{
#' \item{fsalpha}{Alpha parameter for the beta distribution of symptomatic fraction
#' in the likelihood.}
#' \item{fsbeta}{Beta parameter for the beta distribution of symptomatic fraction in the
#' likelihood.  Note the most likely value of fs is \eqn{\frac{\alpha-1}{\alpha+\beta-2}}}
#' }
#'
#' @param prior_weight Factor to multiply the prior by.  Default is 10.  If set
#'   to \code{NULL}, make it equal
#' to the number of counties for which we have confirmed COVID-19 cases.
#' @param fixed_parms A named vector of parameters to use as fixed values for parameters
#' not passed to the likelihood function.  Note these "fixed" parameters can still
#' be overridden by passing an explicit value.
#' @param maxdate Latest date to use in the calibration.
#' @param hparms Named vector of hyperparameters for the likelihood.  See details
#' for supported hyperparameters.
#' @param verbose If \code{TRUE}, have the prior and likelihood print additional
#' diagnostic information.
#' @export
gen_post <- function(prior_weight=10, fixed_parms=NULL, maxdate = NULL, hparms=list(), verbose=FALSE)
{

  hparms <- fill_defaults(as.list(hparms), default_hparms)

  lprior <- gen_prior(hparms, verbose=verbose)
  llik <- gen_likelihood(hparms, fixed_parms, maxdate=maxdate, verbose=verbose)

  if(is.null(prior_weight)) {
    prior_weight <- sum(!is.na(va_county_first_case$firstDay))
  }

  function(parms) {
    logp <- prior_weight * lprior(parms)
    if(verbose) {
      message('prior wgt:\t', prior_weight, '\n----------------\ntotal prior:\t',
              logp, '\n')
    }
    if(is.finite(logp)) {
      logp <- logp + llik(parms)
    }
    logp
  }
}


#' Find the symptomatic infection fraction
#'
#' The fraction of infections that are symptomatic varies over the course of the
#' scenario, so you must pick a range of values to average over.  If no range is
#' specified, the default is days 80-93 (inclusive) from the start of the simulation.
#'
#' @param modrun Model output from \code{\link{run_scenario}}
#' @param trange Range of times to average symptomatic fraction over.  The endpoints
#' are included.
#' @export
fsympto <- function(modrun, trange=c(80,93))
{
  ## silence warnings
  I <- Is <- time <- NULL

  stopifnot(length(trange) == 2)
  if(inherits(trange, 'Date')) {
    trange <- as.integer(trange - as.Date('2020-01-01'))
  }

  ## We have an fs column, but it's computed on a county-by county basis, and
  ## we want the aggregate over the state, which is easiest to get by recomputing
  ## fs from the aggregate data.
  stateagg <-
    dplyr::filter(modrun, time >= trange[1], time <= trange[2]) %>%
    dplyr::group_by(time) %>%
    dplyr::summarise(I=sum(I), Is=sum(Is)) %>%
    dplyr::mutate(fs = Is/(I+Is)) %>%
    dplyr::filter(time == floor(time))         # drop fractional initial time values.

  mean(stateagg$fs)
}

#' Calculate the expected number of hospitalizations for a model
#'
#' Calculate the expected number of hospitalizations by day, based on an assumed
#' fraction of symptomatic patients that eventually go to the hospital.
#'
#' To make this calculation, we have to convert this probability to a rate.  If
#' the recovery rate is \eqn{\gamma = 1/D0}, then the hospitalization rate \eqn{r_h} that gives
#' the desired probability \eqn{P_h} is
#' \deqn{r_h = \gamma \frac{P_h}{1-P_h}}
#'
#' @param ph Per-person hospitalization probability
#' @param modout Model output from \code{\link{run_scenario}}
#' @param mfadjust If \code{TRUE}, adjust for market fraction
#' @return Data frame of time and expected hospitalizations
#' @export
calc_nhosp <- function(ph, modout, mfadjust=TRUE)
{
  rh <- ph/(1-ph) / modout$recoveryTime[1]
  ## We're going to approximate the expected influx of hospital patients as
  ## rh * PopSympto.  Theoretically we should be integrating the hospital population
  ## alongside the other population differential equations, but this approximation is
  ## close enough for our purposes.
  sympt <- modout[c('time','PopSympto')]
  if(mfadjust) {
    sympt[['PopSympto']] <- sympt[['PopSympto']] * modout[['marketFraction']]
  }

  dplyr::group_by(sympt, time) %>%
    dplyr::summarise(expectedHosp = sum(PopSympto * rh))
}

#' Calculate the likelihood adjustment for number of hospital observations
#'
#' The adjustment is
#' \deqn{L_h = \log(\int dp_h Beta(\alpha_h, \beta_h) \prod Pois(N_h, \lambda(p_h)))}
#'
#' @param alpha First shape parameter in the prior for hospitalization fraction
#' @param beta Second shape parameter in the prior for hospitalization fraction
#' @param modout Model output from \code{\link{run_scenario}}
#' @keywords internal
nhosp_likelihood <- function(alpha, beta, modout)
{
  obs <- uva_covid_count[c('time','Admits')]
  modout <- dplyr::filter(modout, time %in% obs$time)

  llsingle <- function(ph) {
    expectHosp <- calc_nhosp(ph, modout)
    cmpHosp <- dplyr::left_join(expectHosp, obs, by='time')
    dbeta(ph, alpha, beta) * prod(dpois(cmpHosp$Admits, cmpHosp$expectedHosp))
  }
  llvector <- function(phvec) {
    sapply(phvec, llsingle)
  }


  log(integrate(llvector, 0,1)$value)
}


#' Generate a likelihood function for single-county simulated data
#'
#' This likelihood function is greatly simplified from the one produced by
#' \code{\link{gen_likelihood}}.  It contains none of the adjustments for
#' hospitalizaton fraction (since simulated data has no hospitalization to
#' compare to), and it doesn't (yet) distinguish between symptomatic and
#' asymptomatic cases.
#'
#' The parameters to this model are beta, A0, D0, Ts, I0, b0, b1, mask_effect,
#' and import_cases.  Mask effect can be omitted; if so, it is fixed at
#' zero.
#'
#' The scenario runs are assumed to start at the earliest time in the observed
#' dataset, even if the time range for comparing to data starts later.  The
#' population is initialized to \code{I0} exposures at that time with the rest
#' of the population initialized as susceptible.
#'
#' Another change in this function is that we use the binomial distribution
#' instead of the hypergeometric.  This saves us some mucking about with
#' ensuring that the susceptible population is greater than the number of
#' positive observations.
#'
#' @param simobs Simulated observations, as returned by \code{\link{simobs}}.
#' @param timerange Length-2 vector giving the earliest and latest dates to use
#'   in the calculation.  If not specified, use all times in the dataset.
#' @param verbose If \code{TRUE}, output additional diagnostic information to
#'   console.
#' @export
gen_simobs_likelihood <- function(simobs, timerange, verbose=FALSE)
{
  t0 <- min(simobs$time)
  pop <- simobs$population[1]

  if(!missing(timerange)) {
    stopifnot(is.numeric(timerange))
    stopifnot(length(timerange == 2))

    simobs <- simobs[simobs$time >= timerange[1] & simobs$time <= timerange[2], ]
  }

  timevals <- simobs$time
  if(timevals[1] != t0) {
    timevals <- c(t0, timevals)
  }

  ## Return this function
  function(p) {
    if(! 'mask_effect' %in% names(p)) {
      p['mask_effect'] <- 0
    }
    stopifnot(setequal(names(p), c('beta','A0','D0','Ts','I0','b0','b1',
                                   'mask_effect', 'import_cases')))
    S0 <- pop-p[['I0']]
    istate <- c(t=t0, S=S0, E=p[['I0']], I=0, Is=0, R=0)
    runout <- run_parmset(p, istate, timevals)
    prevalence <- average_weekly_prevalence(runout, simobs)

    cmp <- dplyr::left_join(prevalence, simobs, by=c('time'))
    if (!all(!is.na(cmp$fi) & !is.na(cmp$ntest))) {
      stop('simobs_likelihood:  Bad values in model output.')
    }
    ## b value declines with more tests, but only to a minimum of 1
    b <- pmax(p[['b0']] - p[['b1']] * log(cmp$ntest), 1)
    cmp$biased_fi <- pmax(padjust(cmp$fi, b), 0)

    x <- cmp$npos
    k <- cmp$ntest
    logl <- dbinom(x, k, cmp$biased_fi, log=TRUE)
    if(any(is.nan(logl))) {
      stop('simobs_likelihood:  Bad output from dbinom')
    }
    sum(logl)
  }
}

#' Generate log-posterior function for filter model
#'
#' @param simobs Simulated observations (see
#'   \code{\link{gen_simobs_likelihood}})
#' @param timerange Time range to use in comparing to observations
#' @param verbose Flag to turn on verbose mode.
#' @export
gen_simobs_posterior <- function(simobs, timerange, verbose=FALSE)
{
  ## default prior function
  logpfun <- default_simobs_prior

  loglfun <-
    if(missing(timerange)) {
      gen_simobs_likelihood(simobs, verbose=verbose)
    }
  else {
    gen_simobs_likelihood(simobs, timerange, verbose)
  }

  function(p) {
    logp <- logpfun(p)
    if(is.finite(logp)) {
      logp + loglfun(p)
    }
    else {
      -Inf
    }
  }
}

### Default prior for parameters in single-locality models.
default_simobs_prior <- function(p) {
  logp <- c(
      dgamma(p['beta'], 10, 30, log=TRUE),
      dgamma(p['A0'], 8, 2, log=TRUE),
      dgamma(p['D0'], 8, 2, log=TRUE),
      dgamma(p['Ts'], 8, 2, log=TRUE),
      dlnorm(p['b0'], 4, 0.2, log=TRUE),
      dlnorm(p['b1'], 0, 1, log=TRUE),
      dlnorm(-p['mask_effect'], -1, 1, log=TRUE),
      dexp(p['import_cases'], 1/50, log=TRUE),
      dexp(p['I0'], 1/50, log=TRUE)
    )
  sum(logp, na.rm=TRUE)
}
