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
      nt <- pmax(1, cmp$ntest)                         # Avoid bad values when ntest == 0
      b <- pmax(oldparms[['b0']] - oldparms[['b1']] * log(nt), 1)
      cmp$biased_fi <- padjust(cmp$fi, b)

      x <- cmp$npos
      m <- ceiling(cmp$biased_fi * cmp$population)
      n <- ceiling((1-cmp$biased_fi) * cmp$population)
      k <- cmp$nt

      dhpenalty <- -1000 * sum(x > m)
      x <- pmin(x, m)

      ## Hypergeometric likelihood.  Ignore if ntest == 0.
      logl <- ifelse(cmp$ntest > 0, dhyper(x, m, n, k, log=TRUE),0)
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
#' @param localityscen \code{LocalScen} object returned from
#'   \code{\link{localize_scenario}}.  If \code{NULL}, then no changes to
#'   parameters will be applied in the future scenario.
#' @return Matrix with a column for time and a column for each state variable.
#' The first row of the matrix (which would contain the initial time) is dropped
#' so that rbinding the output of several consecutive calls to this function will
#' not produce any duplicate rows.
#' @export
run_parmset <- function(parms, istate, timevals, localityscen=NULL)
{
  ode_parms <- c(alpha=1/parms[['A0']], beta=parms[['beta']],
                 gamma=1/parms[['D0']], epsilon=1/parms[['Ts']],
                 mask_effect=parms[['mask_effect']],
                 import_cases=parms[['import_cases']])

  istate <- istate[-c(1)]       # First column is time
  if(is.null(localityscen)) {
    ## Run the ODE solver and drop the initial state (first row)
    rslt <- deSolve::ode(istate, timevals, seir_equations, ode_parms)
    rslt[-c(1), ,drop=FALSE]
  }
  else {
    ## We are going to have parameters changing along the way, so we will need
    ## to run in steps.
    tstop <- max(timevals)              # final stop time
    timebreaks <- scenario_change_times(localityscen)
    timebreaks <- timebreaks[timebreaks < tstop]
    parm_sched <- scenario_parm_value(localityscen,
                                      c(min(timevals), timebreaks),
                                      ode_parms)
    if(max(timebreaks) < tstop) {
      timebreaks <- c(timebreaks, tstop)
    }

    runs <- list()
    length(runs) <- length(timebreaks)
    for(i in seq_along(timebreaks)) {
      ## Find ODE output times for this segment, and remove them from the list
      ## of all output times
      tfin <- timebreaks[i]
      tv <- timevals[timevals <= tfin]
      timevals <- timevals[timevals >= tfin]
      if(max(tv) < tfin) {
        ## Make sure we run right up to the change time, even if it isn't one of
        ## the times requested
        tv <- c(tv, tfin)
        dropfinal <- TRUE
      }
      else {
        dropfinal <- FALSE
      }

      ## Set up parameters for this segment
      for(parm in colnames(parm_sched)) {
        ode_parms[[parm]] <- parm_sched[i, parm]
      }

      ## Run ODE solver
      rr <- deSolve::ode(istate, tv, seir_equations, ode_parms)

      ## Drop first row (it's a repeat of the initial state)
      rr <- rr[-c(1), , drop=FALSE]

      ## Get the new initial state (last row), and drop the last row if it
      ## wasn't wanted.
      istate <- rr[nrow(rr), -c(1)]     # drop the time column.
      if(dropfinal) {
        rr <- rr[-c(nrow(rr)), ]
      }

      runs[[i]] <- rr
    }
    do.call(rbind, runs)
  }
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

#' Run compartment model from an initial state
#'
#' This is a wrapper for \code{\link{run_parmset}} for cases where we have a set
#' of parameters including an initial infection level, and we want to run it
#' from the start of a scenario to a specified time.
#'
#' @param parms Named vector of model parameters: beta, A0, D0, Ts, I0, and
#'   mask_effect.
#' @param obsdata Data frame of observed data.  The only columns used in this
#'   case are time and total population.
#' @param t1 End time for the scenario.
#' @param dt Time step for scenario output.  Default is dt=1.  Note that if
#'   \code{dt > 1}, the end time \code{t1} might not be included in the output.
#' @export
initialize_parmset <- function(parms, obsdata, t1, dt=1)
{
  t0 <- min(obsdata$time)                    # Start at the beginning of the observations
  totpop <- obsdata$population[1]            # Assumed to be constant over time

  timevals <- seq(t0, t1, dt)
  E0 <- parms[['I0']]
  S0 <- totpop - E0
  istate <- c(time=t0, S=S0, E=E0, I=0, Is=0, R=0)
  run_parmset(parms, istate, timevals)
}


### Default prior for filter models
bayes_filter_default_prior <- function(beta0, beta1, imp0, imp1, dt)
{
  ## each week the change in beta is <= .1, 95% of the time (normally distributed)
  wk <- dt/7
  betasig1wk <- 0.05
  betasig <- betasig1wk * sqrt(wk)

  ## Not sure what the best prior is for importation.
  ## For now we'll just go with an exponential distribution with an average of
  ## 2.5/day.
  ## We won't treat multiple weeks as a sum (which would give us a very narrow
  ## result for a multiple-week jump).  Instead, we just multiply the scale
  ## length by the number of weeks.
  imprate <- 1/2.5 * wk


  as.numeric(dnorm(beta1, beta0, betasig, log=TRUE) + dexp(imp1, imprate, log=TRUE) +
             dgamma(beta1, 5, 15, log=TRUE))
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
                       history = TRUE, ids = NULL)
{
  ## Filter to just the times requested
  time <- obsdata[['time']]
  obsdata <- obsdata[time >= tinit & time <= tfinal, ]
  time <- obsdata[['time']]

  ## Vectors for start and end times for filter steps
  strttimes <- time[time < max(time)]
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

  ## Convert to df and add columns for prevalence and expected observations
  addprev <- function(r) {
    r <- as.data.frame(r)

    ## prevalence
    fi <- (r[['I']] + r[['Is']]) / totpop
    r[['fi']] <- fi

    ## expected number of cases
    t <- r[['time']][1] # all the times are the same value, so we just need one.
    iobs <- match(t, obsdata[['time']])
    ntest <- pmax(1, obsdata[['ntest']][iobs])
    b <- r[['b0']] - r[['b1']]*log(ntest)
    ro_biased <- b * fi/(1-fi)
    r[['ncase']] <- ntest * ro_biased/(1+ro_biased)

    ## R0 and Rt values
    r0 <- r[['beta']] * r[['D0']]
    rt <- r0 * r[['S']] / totpop
    r[['R0']] <- r0
    r[['Rt']] <- rt

    r
  }
  rslt <- lapply(rslt, addprev)

  ## If detailed history was requested, add it now.
  if(!is.null(history)) {
    histvars <- c('time', 'S', 'E', 'I','Is','Itot', 'R',
                  'fi', 'ncase', 'R0', 'Rt',
                  'beta','import_cases', 'id')
    if(is.data.frame(history)) {
      stopifnot(setequal(names(history), histvars))
    }
    else {
      history <- tibble::tibble(time=numeric(0), S=numeric(0), E=numeric(0),
                                I=numeric(0), Is=numeric(0), Itot=numeric(0),
                                R=numeric(0), ncase=numeric(0),
                                R0=numeric(0), Rt=numeric(0),
                                fi=numeric(0), beta=numeric(0),
                                import_cases=numeric(0), id=integer(0))
    }
    if(is.null(ids)) {
      ids <- seq(1,nrow(initstates))
    }
    history <- dplyr::bind_rows(history,
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
                 modelfit=fit_table, history=history, ids=ids, obsdata=obsdata),
            class = c('filter-fit','list'))
}


#' Continue a run from a fitted filter model.
#'
#' Unlike \code{\link{fit_filter}}, this function just continues running the
#' model forward in time, with no parameter updates.  Final state and summary
#' statistics are computed and appended to the ones computed for the original
#' fit.  History is appended if it is present in the original fit.
#'
#' @param filterfit Structure of class filter-fit returned by
#'   \code{\link{fit_filter}}.
#' @param tfinal Final time for the projection.
#' @param dt Time step in days.  Note that if \code{dt > 1}, then \code{tfinal}
#'   will not necessarily be in the final output.
#' @param tvintage1 Time for the first projection vintage.  Having multiple
#'   vintages is only possible if history information was included in
#'   \code{filterfit}.
#' @param dtvintage Time step between vintages.
#' @return Data frame with projections in it.
#' @export
project_filter_model <- function(filterfit, tfinal, dt=1, tvintage1=NA, dtvintage=7)
{
  stopifnot(inherits(filterfit, 'filter-fit'))
  tlast <- filterfit[['time']]          # Last time of the model fit.
  if(is.na(tvintage1)) {
    tvintage1 <- tlast
  }
  stopifnot(tvintage1 == tlast || (is.data.frame(filterfit[['history']]) &&
                                   tvintage1 <= max(filterfit[['history']][['time']])))

  vintages <- seq(tvintage1, tlast, dtvintage)


  finalstate <- filterfit[['finalstate']]
  if(is.data.frame(filterfit[['history']])) {
    allstates <- filterfit[['history']]
  }
  else {
    allstates <- as.data.frame(finalstate)
    allstates[['id']] <- seq(1,nrow(allstates))
  }

  totpop <- sum(finalstate[1,c('S','E','I','Is','R')])

  statevars <- c('S','E','I','Is','R')
  parmvars <- c('beta','import_cases')
  parmconst <- c('D0','A0','Ts','mask_effect')

  vintage_runs <-
    foreach::foreach(tinit=vintages, .combine=rbind) %do% {
      timevals <- seq(tinit, tfinal, dt)

      istates <- allstates[allstates[['time']] == tinit, ]

      foreach::foreach(iid=istates[['id']], .combine=rbind) %dopar% {
        irow <- which(istates[['id']] == iid)
        istate <- as.matrix(istates[irow, ])
        vars <- c(tinit, istate[1, statevars])
        parms <- c(istate[1, parmvars], finalstate[irow, parmconst])
        runout <- run_parmset(parms, vars, timevals)
        vintage <- rep(tinit, nrow(runout))
        id <- rep(iid, nrow(runout))
        cbind(id, runout, vintage)
      }
    }
  ## convert to df and summarize.
  vintage_tbl <- as.data.frame(vintage_runs)
  vintage_tbl[['Itot']] <- vintage_tbl[['I']] + vintage_tbl[['Is']]
  vintage_tbl[['fi']] <- vintage_tbl[['Itot']] / totpop

  rslt <-
  dplyr::group_by(vintage_tbl, vintage, time) %>%
  dplyr::group_modify(function(df, grp) {
    as.data.frame(rbind(c(statcols(df, 'S'), statcols(df, 'I'),
                          statcols(df, 'Is'), statcols(df, 'Itot'),
                          statcols(df, 'fi'))))
  })
  attr(rslt, 'locality') <- filterfit[['obsdata']][['locality']][1]
  rslt
}


#### Helper functions used in fit_filter
### Produce median and confidence interval values for a given variable and name
### appropriately.
statcols <- function(tbl, colname)
{
  probs=pnorm(c(-1,0,1))
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
                    R0stats <- statcols(rslt, 'R0')
                    Rtstats <- statcols(rslt, 'Rt')
                    row <- c(time=t, betastats, impstats, fistats, R0stats, Rtstats)
                    if('ncase' %in% names(rslt)) {
                      row <- c(row, statcols(rslt, 'ncase'))
                    }
                    as.data.frame(rbind(row))
                  })
  dplyr::bind_rows(rlist)
}


