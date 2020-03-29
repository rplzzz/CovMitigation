
#' Differential equations for the SIR model
#' 
#' Computes the derivatives of the susceptible (S), infected (I), and recovered (R) 
#' populations.
#' 
#' The variables for this model are S, I, and R, the susceptible, infected, and recovered 
#' populations.  They should be passed as a named vector, in that order.
#' 
#' The parameters for the model are beta, the infection parameter, and gamma, the recovery
#' parameter.  These should be passed in as a named vector; the order doesn't matter.  
#' 
#' Optionally, instead of a number, beta and gamma may be data frames with two columns,
#' 'time' and 'value'.  The time column must start at zero and be strictly
#' increasing, while the beta column may hold any positive values. In this case,
#' beta is considered to be piecewise constant; each time t reaches the the next
#' value of 'time', the value of beta changes to the corresponding value from
#' the beta column.
#' 
#' @param t Simulation time value
#' @param variables Vector of current variable values (see details)
#' @param parameters List or vector of parameter values (see details)
#' @return A list, as described in \code{\link[deSolve]{ode}}.  In this case we provide only
#' the first element of the list, which is a vector of derivative values.
#' @export
sir_equations <- function(t, variables, parameters)
{
  beta <- getparam(t, parameters[['beta']])
  gamma <- getparam(t, parameters[['gamma']])
  with(as.list(variables), {
    dS <- -beta * I * S
    dI <-  beta * I * S - gamma * I
    dR <-  gamma * I
    return(list(c(dS, dI, dR)))
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
    D0                = 7,       # base infection duration
    beta_schedule = data.frame(time=0, value=1), # schedule for relative changes in infection rate
    duration_schedule = data.frame(time=0, value=1), # schedule for relative changes in infection duration
    
    ## Initial state parameters
    S0                = 1,    # Fraction initially susceptible (the rest are uninfected but immune)
    I0                = 1,    # Initial number of infected
    
    ## Disease progression parameters
    symptoFraction = 0.43,   
    hospFraction  = 0.03,   # MMWR(20.7-31.4%)  range I think likely is half that 10.4 - 15.7% with a mean of 13%
                            # but even this is likely an overestimate so I'll start with half of that 6.5%
                            # It looks like we halved this again.  Consult JV.  -rpl
    ratioICUtoAcute = 0.238,          # from the same MMWR 121/508 cases admitted went to the ICU- that is 23.8%  of all admissions
                                      # 26.4% from the 201 person case series Risk Factors Associated with ARDS, C Wu, JAMA 3/13/20
    fractionMV      = 0.328,          # from the 201 person case series Risk Factors Associated with ARDS, C Wu, JAMA 3/13/20
    hospFractionNIPPV   = 0.2605,         # from the 201 person case series Risk Factors Associated with ARDS, C Wu, JAMA 3/13/20
    hospFractionIMV     = 0.0675,         # from the 201 person case series Risk Factors Associated with ARDS, C Wu, JAMA 3/13/20
    
    ## Selection and filtering parameters
    selectMarketShare = 0.2,  # Minimum market share for localities to be in the UVA catchment
    typeAMCDRG = "fractionAllUVA",   # other valid values are "fractionAllVCU"  "fractionRespUVA" "fractionRespVCU" "fractionAllUVA"


    ## time lags (all times in days)
    timeToSympto      = 3,
    timeToAcuteHosp   = 5,
    timeToICU         = 6,
    timeToMV          = 6,
    timeToNIPPV       = 6,
    timeToIMV         = 6
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
#' \item{New community infections in UVA referral base}
#' \item{New symptomatic community infections in UVA referral base}
#' \item{New acute care admissions for UVA}
#' \item{New ICU admissions for UVA}
#' \item{New invasive mechanical ventilation admissions for UVA}
#' }
#' 
#' The parameters for the model are:
#' \describe{
#' \item{T0}{(scalar) Base doubling time.  Doubling time for number of cases when the number of infections
#' is small compared to the total population.}
#' \item{D0}{(scalar) Base recovery time.  Average base recovery time for an infected person.}
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
run_scenario <- function(tmax, params=list(), scenarioName = 'HospCensus') {
  ## Check the parameters and supply defaults as required.
  validate_params(params)
  params <- complete_params(params)
  
  ## Filter the counties to just the ones that meet the market share requirement
  countySelection <- dplyr::filter(marketFractionFinal, 
                                   TypeAMCDRG == params[['typeAMCDRG']],
                                   marketShare >= params[['selectMarketShare']])
  countySelection$Locality <- as.character(countySelection$Locality)
  
  ## Map the single-county function onto our list of counties
  inpatientEstimates <- 
    dplyr::bind_rows(
      mapply(run_single_county, countySelection$Locality, countySelection$marketShare, 
             MoreArgs = list(tmax=tmax, params=params), SIMPLIFY=FALSE)
    )
  
  ## return value:  sum up over the counties
  dplyr::group_by(inpatientEstimates, time) %>%
    dplyr::summarise(newCases = sum(newCases),
                     symptoInfection = sum(symptoInfection),
                     acuteHosp = sum(acuteHosp),
                     icuHosp = sum(icuHosp),
                     imvHosp = sum(IMVHosp),
                     ## The daysToX values should be the same for all localities.  We could 
                     ## retain these columns by including them in the groupings, 
                     ## but on the off chance that they're different for some reason, we'll 
                     ## average them instead.
                     daysToSympto = mean(daysToSympto),
                     daysToAcuteHosp = mean(daysToAcuteHosp),
                     daysToicuHosp = mean(daysToicuHosp),
                     daysToAnyMV = mean(daysToAnyMV),
                     daysToNIPPVHosp = mean(daysToNIPPVHosp),
                     daysToIMVHosp = mean(daysToIMVHosp)
                     ) %>%
    ## Add some identifiers for the scenario
    dplyr::mutate(scenario=scenarioName, doublingTime=params$T0,
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
  ## Need to get the total population from the vaMyAgeBands dataset
  pop <- vaMyAgeBands$Total[vaMyAgeBands$Locality == locality]
  
  ## Calculate derived parameters
  N <- (pop-params$I0) * params$S0
  beta0 <- (1 + (log(2) * params$D0/params$T0)) / (N*params$D0)
  
  beta_schedule <- params$beta_schedule
  beta_schedule$value <- beta_schedule$value * beta0
  
  gamma_schedule <- params$duration_schedule
  gamma_schedule$value <- 1/(gamma_schedule$value * params$D0)
  
  ode_params <- list(beta=beta_schedule, gamma=gamma_schedule)
  initvals <- c(S=N, I=params$I0, R=(pop-params$I0)*(1-params$S0))
  
  timevals <- seq(0, tmax)
  rslt <- as.data.frame(deSolve::ode(initvals, timevals, sir_equations, ode_params))
  
  ## Add some identifiers for the locality
  rslt$locality <- locality
  rslt$population <- pop
  rslt$marketFraction <- mktshare
  
  ## Calculate new cases as lagged difference in S (can't take the difference in I because it is 
  ## also affected by recoveries)
  newcases <- -diff(rslt$S)
  ## This will have one fewer results than the original vector; deem the new cases at t=0 to be 0
  rslt$newCases <- c(0, newcases)
  
  ## Calculate hospitalization numbers and days to hospitalization etc
  rslt$symptoInfection <- rslt$newCases * params$symptoFraction
  rslt$daysToSympto    <- rslt$time + params$timeToSympto
  rslt$acuteHosp       <- rslt$symptoInfection * params$hospFraction * (1-params$ratioICUtoAcute)
  rslt$daysToAcuteHosp <- rslt$time + params$timeToAcuteHosp
  rslt$icuHosp         <- rslt$symptoInfection * params$hospFraction * params$ratioICUtoAcute
  rslt$daysToicuHosp   <- rslt$time + params$timeToICU
  rslt$anyMV           <- rslt$symptoInfection * params$hospFraction* params$fractionMV
  rslt$daysToAnyMV     <- rslt$time + params$timeToMV
  rslt$NIPPVHosp       <- rslt$symptoInfection *  params$hospFraction * params$hospFractionNIPPV
  rslt$daysToNIPPVHosp <- rslt$time + params$timeToNIPPV
  rslt$IMVHosp         <- rslt$symptoInfection  *  params$hospFraction * params$hospFractionIMV
  rslt$daysToIMVHosp   <- rslt$time + params$timeToIMV
  
  rslt
}
