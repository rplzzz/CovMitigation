#' Compute input for the hospital census model from a filter model
#'
#' Run the ensemble members for all counties and aggregate their results,
#' weighted by market share.
#'
#' If the model history is included, it will be prepended onto the projection.
#' If not, the output will start with the projection.
#'
#' The scenario structure, if present, should be a data frame that describes a
#' schedule of hypothetical parameter changes in the future.  Changes for
#' times in the past (prior to the start of the projection) are ignored;
#' however, the most recent such change will take
#' effect in the first future time step.  The format for the scenario structure
#' is a data frame with the following columns:
#' \describe{
#' \item{locality}{The locality the change applies to.  A \code{NA} value here
#' means the change is a default that will apply to any locality without a
#' specified change.}
#' \item{time}{The time at which a change in a parameter occurs.  The parameter
#' will remain at its changed value until another change is imposed.}
#' \item{parm}{Parameter affected by the change}
#' \item{value}{The change factor.  This will be treated as a multiplier on the
#' parameter's value at the end of the historical period.  For example, a value
#' of 1.2 will produce a 20% increase over the parameter's value in the last
#' historical week.}
#' }
#'
#' @section TODO:
#'
#' * Refactor common code between this and \code{\link{project_filter_model}}.
#' * Add options for scenarios with transmissibility changes.
#'
#' @param fitdir Name of directory with saved filter fit objects
#' @param tfinal Max run time for the projections
#' @param scenario_name Name of the scenario.  Defaults to "BAU" (business as
#'   usual).
#' @param mktmin Market share threshold.  Counties where we have a market share
#'   below the threshold will not be run.
#' @param scenario Scenario data structure returned from
#'   \code{\link{Scenario}}.  If a data frame is passed, it will be converted
#'   with a call to \code{\link{Scenario}}.
#' @param tstrt Start of the projection period.  If omitted, the projection will
#'   start at the end of the model fit.  Otherwise, the projection will start at
#'   the last historical time step prior to the requested start time.  Since the
#'   model fitting is done at weekly intervals, this may not be exactly equal to
#'   \code{tstrt}.
#' @return Data frame suitable for input into census model. All case counts are
#'   market share weighted. Columns:
#' \describe{
#' \item{date}{Simulation date}
#' \item{time}{Simulation date expressed as days since Jan 01}
#' \item{scenario}{Ensemble member ID number}
#' \item{popmkt}{Market share weighted population}
#' \item{newCases}{Daily new cases}
#' \item{newSympto}{Daily new symptomatic cases}
#' \item{PopSympto}{Symptomatic population}
#' \item{PopInfection}{Infected population}
#' \item{fracInfection}{Infection prevalence.  Technically this is market share
#'   weighted, but you get the same answer whether you weight by market share or
#'   not.}
#' \item{PopCumulInfection}{Cumulative number of infections}
#' \item{fracCumulInfection}{Cumulative fraction of population infected.  This
#'   value is the same whether or not you weight by market share}
#' }
#' @export
census_model_output <- function(fitdir, tfinal, scenario_name = 'BAU',
                                mktmin=0.1, scenario=NULL, tstrt=NULL)
{
   if(!is.null(tstrt)) {
     stop('Setting tstrt not yet implemented.')
   }

   mfrac <- marketFractionFinal[marketFractionFinal[['TypeAMCDRG']] ==
                                'fractionAllUVA', ]
   counties <- mfrac[mfrac[['marketShare']] >= mktmin, 'Locality']
   ## split vector into variables and parameters
   varnames <- c('time','S','E','I','Is','R')
   parmnames <- c('beta', 'import_cases', 'D0','A0','Ts','mask_effect')

   run_county <- function(locality) {
     if(!is.null(scenario)) {
       localityscen <- scenario_extract_locality(locality, scenario)
     }
     else {
       localityscen <- NULL
     }

     ## load county data
     filename <- file.path(fitdir,paste0('filter-fit-',locality,'.rds'))
     stopifnot(file.exists(filename))
     filterfit <- readRDS(filename)

     ## TODO:  The rest of this function should be refactored.
     dt <- 1                            # to aid in later refactoring.
     tlast <- filterfit[['time']]
     timevals <- seq(tlast, tfinal, dt)
     finalstates <- filterfit[['finalstate']]
     stopifnot(all(finalstates[,'time'] == tlast))

     ids <- filterfit[['ids']]

     ensruns <- foreach::foreach(iens = seq_along(ids),
                                 .combine=rbind, .inorder=FALSE,
                                 .export=c('varnames', 'parmnames')
                                 ) %dopar% {
       istate <- finalstates[iens, ]
       vars <- c(istate[varnames])
       parms <- istate[parmnames]
       rslt <- as.data.frame(run_parmset(parms, vars, timevals, localityscen))  # columns: time, S, E, I, Is, R

       ## prepend history if available
       if(!is.null(filterfit[['history']])) {
         hist <- filterfit[['history']]
         histid <- hist[['id']]
         hist <- hist[histid == ids[iens], varnames]
         ## FML.  The history only has data every 7 days because that's when we
         ## had observational data to compare to.  For now we'll interpolate to
         ## get the days in between.  Ideally we would rerun the model from the
         ## start.
         tinterp <- seq(min(hist[['time']]), max(hist[['time']]))
         histinterp <- data.frame(time=tinterp)
         varnamesxtime <- varnames[-which(varnames == 'time')]
         for(var in varnamesxtime) {
           histinterp[[var]] <- approx(hist[['time']], hist[[var]], tinterp)[['y']]
         }
         rslt <- dplyr::bind_rows(histinterp, rslt)
       }
       else {
         ## If no history available, prepend the final state so that the new
         ## case values are correct in the first time step.
         hist <- as.data.frame(finalstates[iens, , drop=FALSE])
         rslt <- dplyr::bind_rows(hist, rslt)
       }

       ## Compute the columns we need for the census model.
       ## See notes in the run_single_county function (which works for the
       ## static model) for how the new cases are calculated.
       ## TODO:  refactor new sympto calculation so we're not duplicating coe
       ## between here and run_single_county
       ras <- rslt$I / rslt$Is
       deltaIs <- c(0, diff(rslt$Is))
       deltaR <- c(0, diff(rslt$R))
       deltaRs <- deltaR/(1+ras)
       rslt$newSympto <- deltaIs + deltaRs

       rslt$PopInfection <- rslt$I + rslt$Is
       rslt$PopRecovered <- rslt$R
       rslt$PopSympto <- rslt$Is
       rslt$PopSuscept <- rslt$S

       newcases <- diff(rslt$PopInfection) + diff(rslt$PopRecovered)
       rslt$newCases <- c(0, newcases)
       rslt$PopCumulInfection <- cumsum(rslt$newCases)

       ## Add scenario column
       rslt[['scenario']] <- ids[iens]

       ## return result.  Drop the first row because the newSympto entry is not
       ## accurate there.
       rslt[-c(1), c('time', 'scenario', 'newCases', 'newSympto', 'PopSympto', 'PopInfection',
                     'PopRecovered', 'PopCumulInfection', 'PopSuscept')]
     }

     ## Add population and market share to runs
     irow <- which(mfrac[['Locality']] == locality)
     mktfrac <- mfrac[irow, 'marketShare']

     irow <- which(vdhcovid::valocalities[['locality']] == locality)
     pop <- vdhcovid::valocalities[['population']][irow]

     ensruns[['mktfrac']] <- mktfrac
     ensruns[['popmkt']] <- pop*mktfrac

     ## Adjust the other outputs for market fraction
     for(var in c('newCases', 'newSympto', 'PopSympto', 'PopInfection',
                  'PopRecovered', 'PopCumulInfection', 'PopSuscept')) {
       ensruns[[var]] <- ensruns[[var]] * mktfrac
     }

     ensruns
   }

   ## Run all counties
   allensruns <- dplyr::bind_rows(lapply(counties, run_county))

   ## Aggregate across counties
   aggens <-
       dplyr::group_by(allensruns, time, scenario) %>%
       dplyr::summarise(popmkt=sum(popmkt), newCases=sum(newCases),
                        newSympto=sum(newSympto), PopSympto=sum(PopSympto),
                        PopInfection=sum(PopInfection),
                        PopRecovered=sum(PopRecovered),
                        PopCumulInfection=sum(PopCumulInfection),
                        PopSuscept=sum(PopSuscept))
   ## Add the fractional statistics
   aggens[['fracInfection']] <- aggens[['PopInfection']] / aggens[['popmkt']]
   aggens[['fracCumulInfection']] <- aggens[['PopCumulInfection']] / aggens[['popmkt']]
   aggens[['date']] <- aggens[['time']] + as.Date('2020-01-01')

   aggens
}

