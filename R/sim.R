
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
  gamma <- getparam(t, parameters[['gamma']])
  alpha <- getparam (t, parameters[['alpha']])
  epsilon <- parameters[['epsilon']]
  with(as.list(variables), {
    N <- S + E + I + Is + R
    Itot <- I + Is
    expos <- beta * Itot * S / N
    dS <- -expos
    dE <- expos  - alpha*E
    dI <-  alpha*E - gamma * I - epsilon*I
    dIs <- epsilon*I - gamma*Is
    dR <-  gamma * Itot
    return(list(c(dS, dE, dI, dIs, dR)))
  })
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

param_defaults <- 
  list(
    ## Epidemiological model parameters
    T0                = 7.4,     # initial doubling time - this will be turned into an infection rate
    D0                = 4,       # base infection duration - this will be turned into a recovery rate
    A0                = 3,       # base incubation time - this will be turned into a progression rate
    Ts                = 3,       # average time to symptom onset, once progressed to infectious state
    beta_schedule = data.frame(time=0, value=1), # schedule for relative changes in infection rate
    duration_schedule = data.frame(time=0, value=1), # schedule for relative changes in infection duration
    prog_schedule = data.frame(time=0, value=1),  # schedule for relative changes in progression rate
    
    ## Initial state parameters
    S0                = 1,    # Fraction initially susceptible (the rest are uninfected but immune)
    E0                = 0,    # Initial number of exposed (number, not fraction)
    I0                = 1,    # Initial number of infected (a number, not a fraction)
    
    ## day-zero parameter
    day_zero = NULL,
    
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
  c(params, param_defaults[!(names(param_defaults) %in% names(params))])[names(param_defaults)]
}

#' Check that a parameter vector has no unknown parameters
#' 
#' Check against the names in the list of default parameters.
#' 
#' @param params Named list of parameter values
#' @keywords internal
validate_params <- function(params)
{
  iunknown <- !(names(params) %in% names(param_defaults))
  if(any(iunknown)) {
    unknowns <- paste(names(params)[iunknown], collapse=', ')
    warning('Parameter list contained these unknown parameters: ', unknowns)
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
#' \item{T0}{(scalar) Base doubling time.  Doubling time for number of cases when the number of infections
#' is small compared to the total population.}
#' \item{D0}{(scalar) Base recovery time.  Average base recovery time for an infected person.}
#' \item{A0}{(scalar) Base incubation time.  Average tie for an exposed person to 
#' become infected.}
#' \item{(various time lags)}{Time after infection for various events (RTS for details)}
#' \item{(a whole slew of others)}{see definition of \code{param_defaults} for definitions}
#' 
#' }
#' 
#' @param tmax Maximum time (days) to run to
#' @param params List of parameter values, see details
#' @param scenarioName Name for the scenario; this will be copied into results
#' @return Data frame with results (see details) over time
#' @importFrom dplyr %>%
#' @export
run_scenario <- function(tmax, params=list(), scenarioName = 'communityInfection') {
  ## Check the parameters and supply defaults as required.
  validate_params(params)
  params <- complete_params(params)
  
  ## Filter the counties to just the ones that meet the market share requirement
  countySelection <- dplyr::filter(marketFractionFinal, 
                                   TypeAMCDRG == params[['typeAMCDRG']])
  countySelection$Locality <- as.character(countySelection$Locality)
  
  ## Map the single-county function onto our list of counties
  inpatientEstimates <- 
    dplyr::bind_rows(
      mapply(run_single_county, countySelection$Locality, countySelection$marketShare, 
             MoreArgs = list(tmax=tmax, params=params), SIMPLIFY=FALSE)
    )
  
    ## Add some identifiers for the scenario
    dplyr::mutate(inpatientEstimates,
                  scenario=scenarioName, 
                  doublingTime=params$T0,
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
#' @param tmax End time for the simulation
#' @param params Model parameter list
#' @keywords internal
run_single_county <- function(locality, mktshare, tmax, params)
{
  ##msg <- paste(c('loc:', 'mktshare:'), c(locality, mktshare))
  ##message(paste(msg, collapse='  '))
  ## Need to get the total population from the vaMyAgeBands dataset
  pop <- vaMyAgeBands$Total[vaMyAgeBands$Locality == locality]
  if(length(pop) == 0) {
    stop('Cannot find locality:  ',locality)
  }
  
  if(is.null(params[['day_zero']])) {
    dzlookup <- va_county_first_case
  }
  else {
    dzlookup <- params[['day_zero']]
  }

  dayzero <- dzlookup$firstDay[dzlookup$Locality == locality]
  if(is.na(dayzero)) {
    dayzero <- 1 + max(dzlookup$firstDay, na.rm=TRUE)  # crude approx for regions where it hasn't started yet
  }
  dayzero <- dayzero - min(dzlookup$firstDay, na.rm=TRUE)  # Adjust to be relative to first county infected
  
  
  ## Calculate derived parameters
  N <- (pop-params$I0) * params$S0
  
  gamma_schedule <- params$duration_schedule
  gamma_schedule$value <- 1/(gamma_schedule$value * params$D0)
  gamma0 <- getparam(0, gamma_schedule)
  
  alpha_schedule <- params$prog_schedule
  alpha_schedule$value <- 1/(alpha_schedule$value * params$A0)
  alpha0 <- getparam(0, alpha_schedule)
  
  beta0 <- (1 + (log(2) * params$D0/params$T0) * (alpha0+gamma0)/alpha0) / (params$D0)
  
  beta_schedule <- params$beta_schedule
  beta_schedule$value <- beta_schedule$value * beta0
  
  epsilon <- 1/params$Ts
  
  
  ode_params <- list(beta=beta_schedule, gamma=gamma_schedule, alpha=alpha_schedule,
                     epsilon=epsilon)
  initvals <- c(S=(N-params$E0-params$I0)*params$S0,
                E=params$E0, 
                I=params$I0,
                Is=0,
                R=(pop-params$I0-params$E0)*(1-params$S0))
  
  timevals <- seq(0, tmax)
  rslt <- as.data.frame(deSolve::ode(initvals, timevals, seir_equations, ode_params))
  
  ## Adjust for day zero.  Remove any times past tmax
  rslt$time <- rslt$time + dayzero
  #message('dayzero= ', dayzero)
  rslt <- rslt[rslt$time <= tmax,]
  
  ## Add some identifiers for the locality
  rslt$locality <- locality
  rslt$population <- pop
  rslt$marketFraction <- mktshare
  
  ## Calculate new cases as lagged difference in S (can't take the difference in I because it is 
  ## also affected by recoveries)
  newcases <- -diff(rslt$S)
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
  rslt$newSympto <- deltaIs + deltaRs
  
  rslt$PopInfection <- rslt$I + rslt$Is
  rslt$PopRecovered <- rslt$R
  rslt$PopSympto <- rslt$Is
  rslt$PopSuscept <- rslt$S
  rslt$PopCumulInfection <- cumsum(rslt$newCases)

  rslt
}
