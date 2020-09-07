#### Functions for handling scenario specifications.

#### TODO: document all of these functions in a single man page.

#' Construct a scenario structure
#'
#' A scenario structure is a list of data frames.  Each list element pertains to
#' one model variable.  The data frame for the variable gives the change in the
#' parameter by time and locality.
#'
#' The input data frame should include the following columns:
#'
#' \describe{
#' \item{locality}{The locality the change applies to.  A \code{NA} value here
#' means the change is a default that will apply to any locality without a
#' specified change.}
#' \item{time}{The time at which a change in a parameter occurs.  The parameter
#' will remain at its changed value until another change is imposed.}
#' \item{parm}{Parameter affected by the change}
#' \item{isrelative}{Flag indicating whether the change is a multiplier
#' (\code{TRUE}) or an adder (\code{FALSE}) to the base value.}
#' \item{value}{The change factor.  This will be treated as an adjustment to the
#' parameter's value at the end of the historical period.}
#' }

#' @section Notes:
#' Although the relative flag is given at each time step (due to the limitations
#'   of the data frame structure), having an adjustment
#'   that is relative at some times and absolute at others is not supported.
#'   A mixed vector for the \code{isrelative} flag will result in the adjustment
#'   being treated as absolute (i.e., additive).
#'
#' Changes for times in the past (prior to the start of the projection) are
#' ignored; however, the most recent such change will take
#' effect in the first future time step.  Thus, one can keep using the same
#'   scenario structure week after week as the model gets updated.
#'
#' @param scen Input data frame describing the scenario.  See details.
#' @return A scenario structure.
#' @export
Scenario <- function(scen)
{
  if(is.Scenario(scen)) {
    return(scen)
  }

  stopifnot(is.data.frame(scen))
  ## Columns that will appear in the result.  The parm name is required in the
  ## input, but will not appear in the result
  rsltnames <- c('locality', 'time', 'isrelative', 'value')
  reqdnames <- c(rsltnames, 'parm')
  stopifnot(all(reqdnames %in% names(scen)))
  stopifnot(is.numeric(scen[['time']]))
  stopifnot(is.numeric(scen[['value']]))
  stopifnot(is.logical(scen[['isrelative']]))

  ## Split by parameter and drop any unneeded columns.
  rslt <- lapply(split(scen, scen$parm),
                 function(d) {d[ , rsltnames]})

  class(rslt) <- c('Scenario', class(rslt))

  rslt
}


#' Test whether an object is a scenario
#'
#' This function is a convenience wrapper around \code{inherits}
#'
#' @param obj Object to be tested
#' @export
is.Scenario <- function(obj)
{
  inherits(obj, 'Scenario')
}

#' Specialize a scenario to a locality
#'
#' Find the entries in a scenario structure that are germane to a particular
#' locality.  This saves having to filter by locality every time you go to look
#' up an adjustment value.
#'
#' @param scen A \code{\link{Scenario}} object.
#' @param locality Name of the locality to specialize to.
#' @return A \code{LocalScen} structure.  This object is no longer a scenario
#'   object because it can no longer be used with an arbitrary locality.
localize_scenario <- function(scen, locality)
{
  stopifnot(is.Scenario(scen))
  rslt <- lapply(scen, localize_scenario_parm, loc=locality)
  oldclass <- class(rslt)
  class(rslt) <- c('LocalScen', oldclass[oldclass != 'Scenario'])

  rslt
}
### helper function for localize_scenario.  Do not call from anywhere else.
localize_scenario_parm <- function(tbl, loc)
{
  tloc <- tbl[['locality']]
  lcol <- which(names(tbl) == 'locality')
  if(loc %in% tloc) {
    rslt <- tbl[!is.na(tloc) & tloc == loc, -lcol]
    rslt[order(rslt[['time']]), ]
  }
  else {
    if(!any(is.na(tloc))) {
      NULL                  # This will be interpreted as "no variation in this parameter"
    }
    else {
      rslt <- tbl[is.na(tloc), -lcol]
      rslt[order(rslt[['time']]), ]
    }
  }
}

#' Compute parameter value by time for a local scenario
#'
#' @param locscen A \code{LocalScen} object
#' @param timevals Times to retrieve parameter values for
#' @param parmbase Named vector of base values of the parameters
#' @return Matrix of modified parameter values, each variable is in a column,
#'   each time in a row.
scenario_parm_value <- function(locscen, timevals, parmbase)
{
  parmvals <-
      sapply(names(locscen),
             function(parm) {
               if(is.null(locscen[[parm]])) {
                 return(rep(parmbase[[parm]], length(timevals)))
               }

               sched <- locscen[[parm]]
               irow <- sapply(timevals, function(t) {
                 rows <- sched[['time']] <= t
                 if(any(rows)) {
                   max(which(rows))
                 }
                 else {
                   NA
                 }
               })

               isrelative <- all(sched[['isrelative']])
               if(isrelative) {
                 parmbase[[parm]] * ifelse(is.na(irow), 1, sched[['value']][irow])
               }
               else {
                 parmbase[[parm]] + ifelse(is.na(irow), 0, sched[['value']][irow])
               }
             })
}


#' Find the times at which a scenario's parameters change
#'
#' Find the change times (i.e., the times at which at least one parameter
#' changes) across all parameters.
#'
#' @param localscen A \code{LocalScen} object
#' @return Vector of times at which at least one parameter changes.
scenario_change_times <- function(localscen)
{
  skp <- sapply(localscen, is.null)     # variables not present
  as.numeric(unique(
    sort(
      unlist(
        lapply(localscen[!skp], function(x){x$time})
      ))))
}

