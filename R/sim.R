
#' Differential equations for the SEIR model
#'
#' Computes the derivatives of the susceptible (S), exposed (E), infected (I),
#' symptomatic (Is), and recovered (R)
#' populations.  This version of the model omits the population birth and death rates.
#'
#' The variables for this model are S, E, I, Is, and R, the susceptible, infected
#' asymptomatic, infected symptomatic, and recovered
#' populations.  They should be passed as a named vector, in that order.  (The order
#' is important because the ODE solver ignores names.)
#'
#' The parameters for the model are:
#' \describe{
#'   \item{alpha}{Progression rate parameter (i.e., the rate at which people move from exposed to infected)}
#'   \item{beta}{Infection rate parameter:  relative rate at which infected people infect susceptible people.}
#'   \item{gamma}{Recovery rate parameter:  relative rate at which infected people recover (whether symptomatic or not)}
#'   \item{epsilon}{Symptom rate parameter:  relative rate at which asymptomatic people develop symptoms.}
#' }
#' These should be passed in as a named vector; the order doesn't matter.  Also,
#' note that previous versions absorbed the 1/N factor into beta, but this one
#' does not.
#'
#' Optionally, instead of a number, alpha, beta and/or gamma may be data frames with two columns,
#' 'time' and 'value'.  The time column must start at zero and be strictly
#' increasing, while the beta column may hold any positive values. In this case,
#' the parameter is considered to be piecewise constant; each time t reaches the the next
#' value of 'time', the value of the parameter changes to the corresponding value from
#' the value column.  (Time varying epsilon is not currently supported.)
#'
#' @param t Simulation time value
#' @param variables Vector of current variable values (see details)
#' @param parameters List or vector of parameter values (see details)
#' @return A list, as described in \code{\link[deSolve]{ode}}.  In this case we provide only
#' the first element of the list, which is a vector of derivative values.
#' @export
seir_equations <- function(t, variables, parameters)
{
  beta <- getparam(t, parameters[['beta']])
  maskadjust <- exp(parameters[['mask_effect']] * mask_indicator.default(t))
  beta <- beta * maskadjust
  gamma <- getparam(t, parameters[['gamma']])
  alpha <- getparam (t, parameters[['alpha']])
  epsilon <- parameters[['epsilon']]
  S <- variables['S']; E <- variables['E']; I <- variables['I']; Is <- variables['Is']; R <- variables['R']
  N <- S + E + I + Is + R
  Itot <- I + Is
  expos <- beta * Itot * S / N + parameters[['import_cases']]
  dS <- -expos
  dE <- expos  - alpha*E
  dI <-  alpha*E - gamma * I - epsilon*I
  dIs <- epsilon*I - gamma*Is
  dR <-  gamma * Itot
  return(list(c(dS, dE, dI, dIs, dR)))
}


#' Get a time variable parameter from a table of step-changes
#'
#' Check to see if a parameter was given as a table of changes over time.  If so,
#' get the value of the parameter for the current time.  Otherwise, return the parameter's
#' constant value.
#'
#' @param t Simulation time value
#' @param param A data frame giving the step changes in the parameter over time, or a single
#' constant parameter value
#' @return The current value for the parameter, if time variable, or the constant value, if not.
#' @keywords internal
getparam <- function(t, param)
{
  if(is.data.frame(param)) {
    tmin <- min(param$time)
    stopifnot(t >= tmin)
    irow <- max(which(t >= param$time))
    param$value[irow]
  }
  else {
    param
  }
}

### Constant for the day we start tracking the infection, taken to be day 30 (2020-01-31).
infection_t0 <- 30

param_defaults <-
  list(
    ## Epidemiological model parameters
    D0                = 4,       # contagious period - this will be turned into a recovery rate
    A0                = 3,       # base incubation time - this will be turned into a progression rate
    Ts                = 3,       # average time to symptom onset, once progressed to infectious state
    betaSchedule = data.frame(time=0, value=1), # schedule for relative changes in infection rate
    duration_schedule = data.frame(time=0, value=1), # schedule for relative changes in infection duration
    prog_schedule = data.frame(time=0, value=1),  # schedule for relative changes in progression rate
    mask_effect       = 0,       # log-mask effect on transmission
    import_cases      = 0,       # number of imported cases per day.

    ## Initial state parameters
    S0                = 1,    # Fraction initially susceptible (the rest are uninfected but immune)
    E0                = 0,    # Initial number of exposed (number, not fraction)
    I0                = 1,    # Initial number of infected (a number, not a fraction)

    ## future mobility scenario
    mobility_scenario = NULL,

    ## Selection and filtering parameters
    typeAMCDRG = "fractionAllUVA"   # other valid values are "fractionAllVCU"  "fractionRespUVA" "fractionRespVCU" "fractionAllUVA"
  )

#' Add default parameters for parameters not overridden in the input
#'
#' Any parameters not mentioned in the parameter list will have an entry added using the defaults
#' above.
#'
#' @param params Named list of parameter values
#' @return Named list that gives all parameters an explicit value.
#' @keywords internal
complete_params <- function(params)
{
  c(params, param_defaults[!(names(param_defaults) %in% names(params))])
}

#' Check that a parameter vector has no unknown parameters
#'
#' Check against the names in the list of default parameters.
#'
#' @param params Named list of parameter values
#' @keywords internal
validate_params <- function(params)
{
  isbeta <- extract_beta_parms(names(params), 'logical')
  betas <- names(params)[isbeta]
  commons <- names(params[!isbeta])
  iunknown <- !(commons %in% names(param_defaults))

  if(any(iunknown)) {
    unknowns <- paste(commons[iunknown], collapse=', ')
    warning('Parameter list contained these unknown parameters: ', unknowns)
  }

  beta_counties <- beta_parm_loc_extract(betas)
  iunknowncounties <- !(beta_counties %in% va_county_first_case$Locality)
  if(any(iunknowncounties)) {
    unknown_counties <- paste(beta_counties[iunknowncounties], collapse=', ')
    warning('Parameter list contained betas for these unknown counties: ',
            unknown_counties)
  }

  ## Check that parameters have valid values.  These checks aren't comprehensive, just trying to catch
  ## some of the ones that will likely produce obscure errors.
  if('typeAMCDRG' %in% names(params)) {
    stopifnot(!(params[['typeAMCDRG']] %in%
                  c( "fractionAllUVA", "fractionRespUVA", "fractionAllVCU", "fractionRespVCU"
                  )))
  }
}

#' Run a single hospital census scenario
#'
#' Run a hospital census scenario for a single set of model parameters.  Results for several
#' projected variables will be returned in a single data frame.
#'
#' The results returned will be:
#' \itemize{
#' \item{New community infections in the selected counties}
#' \item{New symptomatic community infections in the selected counties}
#' \item{Total infected population in the selected counties}
#' \item{Symptomatic infected population in the selected counties}
#' \item{Recovered population in the selected counties}
#' \item{Remaining susceptible population in the selected counties}
#' \item{Cumulative infected population (includes those who later recovered)}
#' \item{Cumulative infected share of population (includes recovered)}
#' \item{UVA market share of symptomatic infections}
#' }
#'
#' The parameters for the model are:
#' \describe{
#' \item{eta}{(scalar) base value for initial transmissivity}
#' \item{xi}{(scalar) coefficient for population density in initial transmissivity}
#' \item{D0}{(scalar) Base recovery time.  Average base recovery time for an infected person.}
#' \item{A0}{(scalar) Base incubation time.  Average tie for an exposed person to
#' become infected.}
#' \item{(various time lags)}{Time after infection for various events (RTS for details)}
#' \item{(a whole slew of others)}{see definition of \code{param_defaults} for definitions}
#'
#' }
#'
#' @param timevals  Times to output values for.  If only a single value is provided,
#' it is assumed to be a maximum time (days) to run to, and values will be output
#' in one day steps from 0 to tmax.
#' @param params List of parameter values, see details
#' @param scenarioName Name for the scenario; this will be copied into results
#' @param counties If provided, run just these counties
#' @return Data frame with results (see details) over time
#' @seealso \code{\link{run_parms}}
#' @importFrom dplyr %>%
#' @importFrom foreach %do% %dopar%
#' @export
run_scenario <- function(timevals, params=list(), counties = NULL,
                         scenarioName = 'communityInfection')
{
  ## Check the parameters and supply defaults as required.
  validate_params(params)
  params <- complete_params(params)

  if(length(timevals) == 1) {
    timevals <- seq(0, timevals)
  }

  ## Filter the counties to just the ones that meet the market share requirement
  countySelection <- dplyr::filter(marketFractionFinal,
                                   TypeAMCDRG == params[['typeAMCDRG']])
  countySelection$Locality <- as.character(countySelection$Locality)
  if(is.null(counties)) {
    ## If no counties specified, pick all of the ones for which we have supplied
    ## beta values.
    betas <- extract_beta_parms(names(params))
    counties <- beta_parm_loc_extract(betas)
  }

  ## Map the single-county function onto our list of counties
  inpatientEstimates <-
  foreach::foreach(icounty=seq_along(counties), .combine=dplyr::bind_rows,
                   .inorder=FALSE) %dopar% {
      run_single_county(counties[icounty], countySelection$marketShare[icounty],
                        timevals, params)
  }

  beta <- localbeta(params, inpatientEstimates$locality)
  bdt <- calct0(1/params$A0, beta, 1/params$D0)
  ## Add some identifiers for the scenario
  ## TODO: update these identifiers to include time-dependent effects.
  dplyr::mutate(inpatientEstimates,
                scenario=scenarioName,
                beta=beta,
                popfactor=params$xi,
                baseDoublingTime=bdt,
                recoveryTime=params$D0,
                incubationTime=params$A0,
                symptomTime=params$Ts,
                typeAMCDRG=params$typeAMCDRG) %>%
    dplyr::ungroup()
}

#' Run the hospital census model for a single locality
#'
#' @param locality Name of the locality
#' @param mktshare Market share for the locality
#' @param timevals Vector of output times for the simulation
#' @param params Model parameter list
#' @keywords internal
run_single_county <- function(locality, mktshare, timevals, params)
{
  ## Need to get the total population from the vaMyAgeBands dataset
  #message('Running: ', locality)
  pop <- vaMyAgeBands$Total[vaMyAgeBands$Locality == locality]
  if(length(pop) == 0) {
    stop('Cannot find locality:  ',locality)
  }

  ## Find the start date for this county.  Each county is delayed from the overall
  ## start date by a number of days equal to the difference between its first case
  ## and the first case in the state.
  dayzero <- locality_startdate(locality)

  ## Calculate derived parameters
  N <- (pop-params$I0) * params$S0

  gamma_schedule <- params$duration_schedule
  gamma_schedule$value <- 1/(gamma_schedule$value * params$D0)
  ## the table lookups for these can be expensive, so bypass them if we're not
  ## really using the table
  if(nrow(gamma_schedule) == 1) {
    gamma0 <- getparam(0, gamma_schedule)
    gamma_schedule <- gamma0
  }

  alpha_schedule <- params$prog_schedule
  alpha_schedule$value <- 1/(alpha_schedule$value * params$A0)
  if(nrow(alpha_schedule) == 1) {
    alpha0 <- getparam(0, alpha_schedule)
    alpha_schedule <- alpha0
  }

  ## Ugly text matching on every call.  Devise a better way to look this up.
  beta0 <- localbeta(params, locality)
  #message('\tlocality: ', locality,'\tbeta0= ', beta0)

  betaSchedule <- params$betaSchedule
  betaSchedule$value <- betaSchedule$value * beta0
  if(nrow(betaSchedule) == 1) {
    betaSchedule <- beta0
  }

  mobility_table <- local_mobility(locality, params[['mobility_scenario']])

  epsilon <- 1/params$Ts

  #message('\tbeta= ', betaSchedule)
  ode_params <- list(beta=betaSchedule, gamma=gamma_schedule, alpha=alpha_schedule,
                     epsilon=epsilon, mobility_table=mobility_table, zeta=params$zeta,
                     mask_effect=params$mask_effect, import_cases=params$import_cases)
  initvals <- c(S=(N-params$E0-params$I0)*params$S0,
                E=params$E0,
                I=params$I0,
                Is=0,
                R=(pop-params$I0-params$E0)*(1-params$S0))


  ## Silently drop requested output times that are before the start of the infection
  timevals <- timevals[timevals >= dayzero]
  if(length(timevals) == 0) {
    ## There won't be any output within the range of observed data for this county,
    ## so skip it.
    return(NULL)
  }
  ## Ensure that the calculation starts on the correct date.
  if(min(timevals) > dayzero) {
    timevals <- c(dayzero, timevals)
    dropfirst <- TRUE                # drop the first row from the results, since it wasn't requested.
  }
  else {
    dropfirst <- FALSE
  }

  rslt <- as.data.frame(deSolve::ode(initvals, timevals, seir_equations, ode_params))

  ## Add some identifiers for the locality
  rslt$locality <- locality
  rslt$population <- pop
  rslt$marketFraction <- mktshare

  ## Calculate new cases as lagged difference in total infected pop minus lagged
  ## difference in recovered pop
  newcases <- diff(rslt$I + rslt$Is) + diff(rslt$R)
  ## This will have one fewer results than the original vector; deem the new cases at t=0 to be 0
  rslt$newCases <- c(0, newcases)

  rslt$fs <- rslt$Is / (rslt$Is+rslt$I)

  ## New symptomatic cases are tricky because we have cases entering and leaving
  ## the pool.  Here's how we do it:
  ## deltaIs = newSympto - deltaRs
  ## deltaR  = deltaRs + deltaRa
  ## deltaRa / deltaRs = I/Is = ras
  ##    --> deltaR = (1+ras) deltaRs  --> deltaRs = deltaR/(1+ras)
  ## newSympto = deltaIs + deltaRs
  ras <- rslt$I / rslt$Is
  deltaIs <- c(0, diff(rslt$Is))
  deltaR <- c(0, diff(rslt$R))
  deltaRs <- deltaR/(1+ras)
  rslt$newSympto <- newsympto(rslt$I, rslt$Is, rslt$R)

  rslt$PopInfection <- rslt$I + rslt$Is
  rslt$PopRecovered <- rslt$R
  rslt$PopSympto <- rslt$Is
  rslt$PopSuscept <- rslt$S
  rslt$PopCumulInfection <- cumsum(rslt$newCases)

  ## Drop the initial time, if it was not amongst the times requested
  if(dropfirst) {
    rslt <- rslt[-c(1),]
  }

  rslt
}

#' Run the scenario for a vector of likelihood parameters
#'
#' This function extracts the parameters that need to be passed to
#' \code{\link{run_scenario}}, figures out the time values to use,
#' and filters out the fractional days from the output.
#'
#' @param parms Vector of parameters used in the likelihood function (q.v.,
#' \code{\link{gen_likelihood}}).  Unlike with the likelihood function, there
#' is no provision for default parameters
#' @param scenario_name Optional name for the scenario
#' @param tmax Maximum time to run to (default is 273, which is 2020-09-30)
#' @param aggregate If \code{TRUE}, aggregate across counties at each time.
#' Otherwise, return the unaggregated data.
#' @seealso \code{\link{run_scenario}}
#' @export
run_parms <- function(parms, scenario_name='Scen', tmax=273, aggregate=FALSE)
{
  seirparms <- parms[c('eta', 'xi', 'zeta', 'D0', 'A0', 'I0', 'Ts',
                       'mask_effect', 'mobility_scenario')]

  stopifnot(tmax > infection_t0)
  tvals <- seq(infection_t0, ceiling(tmax))
  modout <- run_scenario(tvals, as.list(seirparms), scenarioName = scenario_name)

  modout <- modout[modout$time - floor(modout$time) < 1e-6,]
  modout
}


#' Calculate incidence of symptomatic cases
#'
#' Given a time series of symptomatic and asymptomatic infections, and recovered
#' population, calculate the new symptomatic cases at each time period.
#'
#' New symptomatic cases are tricky because we have cases entering and leaving
#' the pool.  Basically, we look at the change in Is and back out the portion of
#' the change in R that is attributable to recovery from symptomatic cases.
#' Because we assume the recovery rates for symptomatic and asymptomatic cases
#' are equal, this is proportional to the ratio of symptomatic to asymptomatic
#' infections.
#'
#' Formally:\cr
#' ras = I / Is\cr
#' deltaIs = newSympto - deltaRs\cr
#' deltaR  = deltaRs + deltaRa\cr
#' deltaRa / deltaRs = I/Is = ras\cr
#'    --> deltaR = (1+ras) deltaRs  --> deltaRs = deltaR/(1+ras)\cr
#' newSympto = deltaIs + deltaRs\cr
#' where deltaX is the change in X from the previous time step.
#'
#' Because the results depend on differences from the previous time step, we
#' can't calculate the new cases for the first time step.  If the first time
#' step is the beginning of the infection, then we're looking at an initial
#' steady state, and the initial difference values are therefore zero, so that's
#' what we use.  If the first data point is in the middle of a run, then this
#' assumption is no good, and we have to discard the first value.  The
#' \code{dropfirst} argument allows this to happen automatically.
#'
#' @param I Vector of asymptomatic infections over time.
#' @param Is Vector of symptomatic infections over time.
#' @param R Vector of recovered/removed population over time.
#' @param dropfirst If \code{TRUE}, drop the first entry (see details)
newsympto <- function(I, Is, R, dropfirst=FALSE)
{
  ras <- I/Is
  deltaIs <- c(0, diff(Is))
  deltaR <- c(0, diff(R))
  deltaRs <- deltaR/(1+ras)
  newsympto <- deltaIs + deltaRs
  if(dropfirst) {
    newsympto[-1]
  }
  else {
    newsympto
  }
}
