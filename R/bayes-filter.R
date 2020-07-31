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
#' @param initstates Matrix of initial states and parameter samples (see details).
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
  statevars <- c('time', 'S','E','I','Is','R')
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
                 mask_effect=parms[['mask_effect']],
                 import_cases=parms[['import_cases']])

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


#' Fit a Bayesian filter over time
#'
#' This function calls \code{\link{bayesian_filter}} repeatedly to fit an
#' evolving set of parameters to a dataset.  The result is organized into a
#' table of model state and parameters over time.
#'
#' The \code{history} argument allows us to optionally keep detailed histories
#' for the individual ensemble members, including beta, imports, I, Is, and
#' Itot.  Otherwise, the only history you get is is summary statistics (median
#' and 95% confidence bounds) for those quantities.  To get the full history,
#' either pass in the history output from a previous run (the new history will
#' be appended), or if earlier history is not available, set to any non-null
#' value (the new history will be returned).
#'
#' If history is being kept, then we need to assign serial numbers to the
#' ensemble members.  By default, these are assigned sequentially, but there are
#' reasons why they might need to be specified explicitly.  One is if the
#' ensemble matrix is reordered for some reason.  The other is if some ensemble
#' members are dropped or replaced.
#'
#' @param initstates Matrix containing ensemble of initial states and parameter
#'   values (see description in \code{\link{bayesian_filter}}.
#' @param tinit Initial time corresponding to the states provided in
#'   \code{initstates}
#' @param tfin Final time (i.e., the time of the last observation to use in
#'   fitting the filter)
#' @param obsdata Observed data to use in the fit
#' @param history Detailed history for the individual ensemble members (see
#'   details).
#' @param ids List of serial numbers for the ensemble members.  This is only
#'   used if we are keeping history.  See details for why you might want to
#'   specify this explicitly.
#' @export
fit_filter <- function(initstates, obsdata, tinit, tfinal,
                       history = NULL, ids = NULL)
{
  ## Filter to just the times requested
  time <- obsdata[['time']]
  obsdata <- obsdata[time >= tinit & time <= tfinal, ]
  time <- obsdata[['time']]

  ## Vectors for start and end times for filter steps
  strttimes <- time[time < tfinal]
  endtimes <- time[time > tinit]
  ntime <- length(strttimes)
  stopifnot(length(endtimes) == ntime)

  ## loop over times.  This loop is necessarily sequential
  rslt <- list()
  length(rslt) <- ntime
  locality <- obsdata[['locality']][[1]]
  for(i in seq_along(strttimes)) {
    if(i == 1) {
      ist <- initstates
    }
    else {
      ist <- rslt[[i-1]]
    }

    rslt[[i]] <- bayesian_filter(ist, locality, strttimes[i], endtimes[i],
                                 obsdata)
  }
  finalstate <- rslt[[length(rslt)]]

  ## Figure out the total population
  timecol <- which(colnames(initstates) == 'time')
  totpop <- sum(initstates[1,-c(timecol)])

  ## Convert to df and add a column for prevalence
  addprev <- function(r) {
    r <- as.data.frame(r)
    r[['fi']] <- (r[['I']] + r[['Is']]) / totpop
    r
  }
  rslt <- lapply(rslt, addprev)

  ## If detailed history was requested, add it now.
  if(!is.null(history)) {
    histvars <- c('I','Is','Itot','fi','beta','import_cases', 'id')
    if(is.data.frame(history)) {
      stopifnot(setequal(names(history), histvars))
    }
    else {
      history <- tibble::tibble(I=numeric(0), Is=numeric(0), Itot=numeric(0),
                                fi=numeric(0), beta=numeric(0),
                                import_cases=numeric(0), id=integer(0))
    }
    if(is.null(ids)) {
      ids <- seq(1,nrow(initstates))
    }
    dplyr::bind_rows(history,
                     lapply(rslt, function(r) {
                       ## fi was added, but we didn't keep Itot
                       r$Itot <- r$I + r$Is
                       r$id=as.integer(ids)
                       r[ , histvars]
                     }))
  }

  ## Gather results into a table.  The times for the results are the end times
  ## of each run.
  fit_table <- collate_results(rslt, endtimes)

  ## Return the result, including both final state and the table of fit values
  structure(list(finalstate=finalstate, time=endtimes[length(endtimes)],
                 modelfit=fit_table, history=history),
            class = c('filter-fit','list'))
}


#### Helper functions used in fit_filter
### Produce median and confidence interval values for a given variable and name
### appropriately.
statcols <- function(tbl, colname)
{
  probs=c(0.025, 0.5, 0.975)
  stats <- quantile(tbl[[colname]], probs, names=FALSE)
  names(stats) <- paste0(colname, c('lo','','hi'))
  stats
}

### Collate a list of results into a single data frame.
collate_results <- function(rsltlist, timevals)
{
  rlist <- lapply(seq_along(rsltlist),
                  function(i) {
                    t <- timevals[i]
                    rslt <- rsltlist[[i]]
                    probs <- c(0.025, 0.5, 0.975)
                    betastats <- statcols(rslt, 'beta')
                    impstats <- statcols(rslt, 'import_cases')
                    fistats <- statcols(rslt, 'fi')
                    row <- c(time=t, betastats, impstats, fistats)
                    m <- matrix(row, nrow=1)
                    d <- as.data.frame(m)
                    colnames(d) <- names(row)
                    d
                  })
  dplyr::bind_rows(rlist)
}
