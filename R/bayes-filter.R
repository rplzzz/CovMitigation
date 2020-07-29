#' Project infection in a single county using a Bayesian filter
#'
#' This function runs a single county one week forward in time from a distribution
#' of model states.  For each starting state it produces one sample from the posterior PDF
#' of beta values and importation values for the next week.  The result is a distribution
#' of one-week beta values and importation values, along with a new distribution of
#' states at the end of the week.
#'
#' The initial states should be a matrix with samples in rows and variables in
#' columns.  The variables should be the state variables of the infection equations
#' (S, E, I, Is, and R), plus the beta and import_cases parameters.
#'
#' The observed data should be the first element of the list returned by
#' \code{\link{get_obsdata}}, filtered to the locality being analyzed.  If
#' omitted, we will fetch the data and filter it; however, if this function will
#' be called consecutively for several weeks, running, it will save some time to
#' do the filtering once and pass the table in.
#'
#' @param initstates Vector of initial state and parameter samples (see details).
#' @param tinit Time corresponding to the initial states.
#' @param tfin End time for the forecast.
#' @param prior Prior PDF function.  It should take a vector of five parameters:  start and end
#' beta values, start and end importation values, and size of time step.
#' @param nsamp Number of Monte Carlo Samples to run for each state.  The final set
#' of parameters from the Monte Carlo will be used to generate the finals state.
#' @param obsdata Observed number of tests and number of positive tests by
#'   county. (See details)
#' @return A table of parameter values and final states corresponding to each initial
#' state.  This output is suitable to feed back into the function as \code{initstates}
#' for a subsequent week.
#' @export
bayesian_filter <- function(initstates, locality, tinit, tfin,
                            obsdata=NULL,
                            nsamp = 1000,
                            prior=bayes_filter_default_prior)
{
  ## values independent of the parameters or the state
  dt <- tfin - tinit
  timevals <- seq(tinit, tfin)
  sampscale <- c(0.1, 20)
  statevars <- c('S','E','I','Is','R')
  parmvars <- c('beta', 'import_cases', 'D0', 'A0', 'Ts',
                'mask_effect', 'b0', 'b1')  # Note that we don't need I0.

  if(is.null(obsdata)) {
    obsdata <- get_obsdata()[[1]]
    obsdata <- obsdata[obsdata$locality == locality, ]
  }
  if(!('population' %in% names(obsdata))) {
    ## Grab the population for this locality
    loc <- vdhcovid::valocalities
    pop <- loc[loc$locality == locality, 'population', drop=TRUE]
    obsdata$population <- pop
  }

  ## Loop over the parameter and state samples.
  loopfn <- function(irow) {
    istate <- initstates[irow, statevars]
    oldparms <- initstates[irow, parmvars]

    priorfun <- function(p) {
      ## p[1] == beta; p[2] == import_cases
      prior(oldparms[1], p[1], oldparms[2], p[2], dt)
    }
    likelihoodfun <- function(p) {
      newparms <- oldparms
      newparms[1:2] <- p
      runout <- run_parmset(newparms, istate, timevals)
      prevalence <- average_weekly_prevalence(runout, obsdata)

      cmp <- dplyr::left_join(prevalence, obsdata, by=c('time'))
      stopifnot(all(!is.na(cmp$fi) & !is.na(cmp$ntest))) # All data must be present
      b <- pmax(oldparms[['b0']] - oldparms[['b1']] * log(cmp$ntest), 1)
      cmp$biased_fi <- padjust(cmp$fi, b)

      x <- cmp$npos
      m <- ceiling(cmp$biased_fi * cmp$population)
      n <- ceiling((1-cmp$biased_fi) * cmp$population)
      k <- cmp$ntest

      dhpenalty <- -1000 * sum(x > m)
      x <- pmin(x, m)

      logl <- dhyper(x, m, n, k, log=TRUE)
      sum(logl)
    }

    postfun <- function(p) {
      pr <- priorfun(p)
      if(is.finite(pr)) {
        pr + likelihoodfun(p)
      }
      else {
        pr
      }
    }

    opt <- optim(oldparms[1:2], postfun, control=list(fnscale=-1))

    ms <- metrosamp::metrosamp(postfun, opt$par, nsamp, 1, sampscale)
    if(ms$accept < 0.1) {
      warning('Low sampling acceptance rate: accept = ', ms$accept,
              '  with scale = ', paste(ms$scale, collapse=', '))
    }
    else if(ms$accept > 0.6) {
      warning('High sampling acceptance rate: accept = ', ms$accept,
              '  with scale = ', paste(ms$scale, collapse=', '))
    }

    ## Find the new state for our final parameters
    newparms <- oldparms
    newparms[1:2] <- ms$plast
    runout <- run_parmset(newparms, istate, timevals)

    ilast <- nrow(runout)

    c(runout[ilast, statevars], newparms)

  }
  foreach::foreach(irow=seq(1, nrow(initstates)), .combine=rbind) %dopar%
      loopfn(irow)
}


#' Run the compartment model for a single set of parameters
#'
#' By contrast to \code{\link{run_parms}}, this function runs just the ODE
#' integration for the SEIR equations.
#'
#'
#' @param parms Named vector of model parameters: beta, A0, D0, Ts, and mask_effect
#' @param istate Named vector of initial values for state variable: S, E, I, Is, and R
#' @param timevals Vector of times to output results.  The values of the state
#' variables will be initialized at the first time in the list.
#' @return Matrix with a column for time and a column for each state variable.
#' The first row of the matrix (which would contain the initial time) is dropped
#' so that rbinding the output of several consecutive calls to this function will
#' not produce any duplicate rows.
#' @export
run_parmset <- function(parms, istate, timevals)
{
  ode_parms <- c(alpha=1/parms[['A0']], beta=parms[['beta']],
                 gamma=1/parms[['D0']], epsilon=1/parms[['Ts']],
                 mask_effect=parms[['mask_effect']])

  istate <- istate[-c(1)]       # First column is time
  ## Run the ODE solver and drop the initial state (first row)
  rslt <- deSolve::ode(istate, timevals, seir_equations, ode_parms)
  rslt[-c(1), ,drop=FALSE]
}

average_weekly_prevalence <- function(runout, obsdata)
{
  timecol <- which(colnames(runout) == 'time')
  weeks <- sapply(runout[ , timecol],
                  function(x) {which.max(obsdata[['time']] >= x)})

  totalpop <- sum(runout[1, -c(timecol)])
  prevalence <- (runout[, 'I'] + runout[, 'Is']) / totalpop

  avgprev <- as.numeric(sapply(split(prevalence, weeks), mean))
  data.frame(time=obsdata[['time']][unique(weeks)], fi=avgprev)
}

bayes_filter_default_prior <- function(beta0, beta1, imp0, imp1, dt)
{
  ## each week the change in beta is <= .1, 95% of the time (normally distributed)
  wk <- dt/7
  betasig1wk <- 0.05
  betasig <- betasig1wk * sqrt(wk)

  ## Not sure what the best prior is for importation.  I was thinking
  ## exponential with a scale length of the previous weekly importation rate
  ## (per week, minimum of 10/wk), but I don't think we want to treat multiple
  ## weeks as a sum (which would give us a very narrow result for a
  ## multiple-week jump).  Instead, we just multiply the scale length by the
  ## number of weeks.
  imprate <- 1/pmax(imp0, 10) * wk


  as.numeric(dnorm(beta1, beta0, betasig, log=TRUE) + dexp(imp1, imprate, log=TRUE))
}
