#' Caclulate model coefficients from parameters
#'
#' We parameterize the model SEIR in terms of initial doubling time, incubation period,
#' and recovery time.  The model equations are given in terms of rates,
#' \eqn{\alpha}, \eqn{\beta}, and \eqn{\gamma}.  Also of interest is the basic
#' reproduction parameter \eqn{R_0}.
#'
#' Recovery rate and progression rate are just the reciprocals of their respective
#' periods:
#' \deqn{\alpha = 1/A_0}
#' \deqn{\gamma = 1/D_0}
#'
#' Beta is calculated from \eqn{\alpha}, \eqn{\gamma}, and \eqn{t_d} as
#' \deqn{\beta = \gamma + \ln 2 \frac{(\alpha+\gamma)}{\alpha t_d}}
#'
#' \eqn{R_0} is calculated from \eqn{\beta} and \eqn{\gamma}.  Note that in the
#' absence of mortality during the exposed period, \eqn{R_0} is independent of
#' \eqn{\alpha}.
#' \deqn{R_0 = \frac{\beta}{\gamma}}
#'
#' @param parms Named vector of model parameters
#' @param alpha SEIR alpha parameter
#' @param beta SEIR beta parameter
#' @param gamma SEIR gamma parameter
#' @name coefs
NULL


#' @describeIn coefs Calculate progression rate \eqn{\alpha} from model parameters
#'
#' @export
calcalpha <- function(parms)
{
  1/parms[['A0']]
}

#' @describeIn coefs Calculate transmissibility \eqn{\beta} from model parameters.
#'
#' This is the approximate version, using the Taylor series expansion.
#'
#' @export
calcbeta_approx <- function(parms)
{
  alpha <- calcalpha(parms)
  gamma <- calcgamma(parms)
  td <- parms[['T0']]
  gamma + log(2) *(alpha+gamma) / (alpha*td)
}

#' @describeIn coefs Calculate recovery rate \eqn{\gamma} from model parameters.
#'
#' @export
calcgamma <- function(parms)
{
  1/parms[['D0']]
}

#' @describeIn coefs Calculate beta, given model input parameters
#'
#' This is the exact version, solving the nonlinear equation for beta.
#'
#' @export
calcbeta <- function(parms)
{
  beta_guess <- calcbeta_approx(parms)
  alpha <- calcalpha(parms)
  gamma <- calcgamma(parms)
  targfun <- function(beta) {
    parms[['T0']] - calct0(alpha, beta, gamma)
  }
  solv <- nleqslv::nleqslv(beta_guess, targfun)
  if(solv$termcd != 1) {
    warning('calcbeta: failed to achieve convergence: msg= ', solv$message)
  }
  solv$x
}


extract_beta_parms <- function(pnames, rtn='names')
{
  isbeta <- grepl('^beta_', pnames)

  switch(rtn,
         logical=isbeta,
         indices = which(isbeta),
         pnames[isbeta])
}

beta_parm_loc_extract <- function(bpnames)
{
  ## parameter names are of the form beta_localityName (e.g., beta_AlbemarleCounty)
  substring(bpnames, 6)
}

#' Calculate the local value of beta from population density
#'
#' Our formula is \eqn{\beta = \exp(\eta + \xi d)}, where \eqn{d} is the
#' standardized population density.
#'
#' This function is vectorized over localities.
#'
#' @export
localbeta <- function(parms, localities)
{
  default_beta_val <- 1.0/3.0      # default for any counties unspecified.
  idx <- match(paste0('beta_', localities), names(parms))

  betavals <- parms[idx]
  betavals[is.na(betavals)] <- default_beta_val
  simplify2array(betavals)
}

#' @describeIn coefs Calculate effective reproduction number \eqn{R_e} from model parameters.
#'
#' @export
calcreff <- function(parms)
{
  calcbeta(parms)/calcgamma(parms)
}

#' @describeIn coefs Calculate the growth rate (lplus) given alpha, beta, gamma coefficients
#'
#' @export
calclplus <- function(alpha, beta, gamma)
{
  apg <- alpha+gamma
  0.5*(-apg + sqrt(apg^2 + 4*alpha*(beta-gamma)))
}

#' @describeIn coefs Calculate the doubling time (t0), given alpha, beta, gamma coefficients
#'
#' @export
calct0 <- function(alpha, beta, gamma)
{
  log(2) / calclplus(alpha, beta, gamma)
}

#' Adjust probabilities by a bias factor
#'
#' Given a vector of probabilities and a bias factor (or vector of bias factors),
#' compute the new probabilities adjusted for the bias factor.
#'
#' The bias factor acts as a multiplier on the odds ratio, so
#' \deqn{\frac{p'}{1-p'} = b \frac{p}{1-p}.}  Thus, for example, if \eqn{p=0.1},
#' the odds ratio is \eqn{r=1/9}.  A bias factor of 10 would produce \eqn{r'=10/9},
#' for an ajusted probability of \eqn{p' \approx 0.53}.
#'
#' @param p Vector of probability values
#' @param b Bias factor
#' @export
padjust <- function(p, b)
{
  biased_odds <-  p/(1-p) * b
  biased_odds / (1 + biased_odds)
}

#' Find a name for a backup file
#'
#' The backup name is formed by first cutting off the extension (everything after
#' the final '.') and adding the tag ```_NNN_```, where N counts up from the
#' number of the last backup.  Then the extension is added back.
#'
#' @param name Name(s) of the files to back up
#' @export
namebackup <- function(name)
{
  extbrk <- regexpr('\\.[^.]*$', name)
  namestem <- ifelse(extbrk > 1, substr(name, 1,extbrk-1), name)
  nameext <- ifelse(extbrk > 1, substring(name, extbrk), '')
  backupbrk <- regexpr('_[0-9]+_$', namestem)
  backupindx <- ifelse(backupbrk > 1,
                       1 + as.integer(substr(namestem,
                                             backupbrk+1,
                                             backupbrk + attr(backupbrk, 'match.length')-2)),
                       0)
  backupstem <- ifelse(backupbrk > 1,
                       substr(namestem, 1, backupbrk-1),
                       namestem
                       )

  ## paste the parts together for the return value
  paste0(backupstem, sprintf('_%03d_', backupindx), nameext)
}

#' Aggregate model output to weekly resolution
#'
#' Compute the sum or average of the daily output from the model, for weeks
#' ending on dates specified in the output.  Aggregation will be performed only
#' for the requested columns; other columns will have their values at the
#' specified dates copied over unmodified.
#'
#' @param weekending End dates for the weeks being calculated.
#' @param df Dataset to operate on
#' @param cols Columns to aggregate
#' @param aggfun Function to aggregate with (default is \code{mean})
#' @param nday Number of days in the "week".  Default is 7
#' @export
wkagg <- function(weekending, df, cols, aggfun=mean, nday=7)
{
  doagg <- function(date) {
    dend <- date
    dstrt <- date - nday
    iend <- which(df$date == dend)

    usedays <- df$date > dstrt & df$date <= dend
    if(sum(usedays) == 0) {
      ## no days in the dataset for this week
      return(NULL)
    }
    dlast <- max(df[['date']][usedays])
    iend <- which(df$date == dlast)

    ## Date must be in the data frame exactly once.
    stopifnot(length(iend) <= 1)

    rslt <- list(date=date)
    for(col in names(df)) {
      if(col %in% cols) {
        vals <- df[[col]][usedays]
        rslt[[col]] <- aggfun(vals[!is.na(vals)])
      }
      else {
        rslt[[col]] <- df[[col]][iend]
      }
    }

    tibble::as_tibble(rslt)
  }

  dplyr::bind_rows(lapply(weekending, doagg))
}

#' Compute transmissibility adjustment for mobility
#'
#' Compute multiplicative adjustment to transmissibility, based on percentage
#' a population mobility measurement
#' The adjustment is \code{exp(zeta * mobility)}
#'
#' The mobility table is processed in \code{\link{local_mobility}}.
#'
#' @param t Days since 2020-01-1
#' @param zeta Model coefficient for influence of mobility on transmissibility
#' @param mobility_table Mobility table for the locality being processed (see details)
#' @export
mobility_adjust <- function(t, zeta, mobility_table)
{

  ttbl <- mobility_table[[1]]
  ival <- which.max(ttbl >= t)     # sentinel value guarantees that there is at least one TRUE in the table.
  mobval <- mobility_table[[2]][ival]

  exp(zeta * mobval)
}

#' Extract the mobility table for a specified locality.
#'
#' The mobility data currently in use is the Google mobility data, stored in the
#' \code{\link{va_mobility_daily}} dataset.
#'
#' The mobility dataset has columns for retail, grocery, parks, transit, work, and
#' home categories.  (The latter is actually inversely related to mobility, since
#' it represents staying at home.)  The default is to use the 'home' column, but
#' this can be changed by setting the option \code{CovMitigation.mobility_column}.
#'
#' Setting the mobility category via an option is manifestly a terrible idea, since
#' the mobility column used is not recorded anywhere in the output.  Once we decide
#' which category best captures the effect that we are after, we will set that
#' column as the default and disable the set-by-option capability
#'
#' @param locality Name of the locality to extract
#' @param scenario Future mobility scenario to use.  If \code{NULL}, then use the
#' base scenario (with mobility constant after the last data point).  Otherwise,
#' the name of the future scenario to use.
#' @return A list with two vectors, \code{t} and \code{mobility}, in that
#' order.
#' @export
local_mobility <- function(locality, scenario=NULL)
{
  mobility_col <- 'work'
  if(is.null(scenario)) {
    mobility_base_table <- va_mobility_weekly
  }
  else {
    ## For now we only have one future scenario, so just use that for any non-null
    ## vaue of scenario
    mobility_base_table <- va_future_mobility_weekly
  }
  mobility_table <-
    mobility_base_table[mobility_base_table$locality == locality ,
                      c('t', mobility_col)]
  names(mobility_table) <- c('t', 'mobility')
  mobility_table <- mobility_table[!is.na(mobility_table[['mobility']]),]

  nr <- nrow(mobility_table)
  if(nr == 0) {
    ## No data in the table, so create a dummy table that will produce no
    ## adjustment.
    mobility_table <- tibble::tibble(t=c(0,1e6), mobility=c(0,0))
  }

  as.list(mobility_table)
}

#' Generate a vector of parameters for the posterior PDF
#'
#' Create a vector of parameters for use with the posterior function.  All
#' counties are initialized with the same beta value, and values can be
#' specified for the statewide parameters.  Any parameters not specified will be
#' given default values.
#'
#' @param localities Vector of localities to use in the calculation.  If
#'   omitted, use all localities for which there is case data.
#' @param beta Value to use for transmissibility.  All localities will get the
#'   same value.
#' @param common_parms Values to use for the statewide parameters.  If omitted,
#'   default values will be used
#' @return Named vector of parameter values.
#' @export
gen_parm_vec <- function(localities=NULL, beta=NULL, common_parms=NULL)
{
  cpvals <- c(D0=4, A0=2, Ts=4, I0=30, mask_effect=-1.5, b0=50, b1=1) # default values
  if(!is.null(common_parms)) {
    beta_specified <- extract_beta_parms(names(common_parms), 'logical')
    if(any(beta_specified)) {
      warning('Beta values found in common_parms.  These will be ignored.')
      common_parms <- common_parms[!beta_specified]
    }
    pspec <- names(common_parms)
    cpvals[pspec] <- common_parms
  }

  if(is.null(localities)) {
    localities <-
      va_county_first_case$Locality[!is.na(va_county_first_case$firstDay)]
  }

  if(is.null(beta)) {
    beta <- 1.0/3.0
  }
  if(length(beta) > 1) {
    warning('Vector supplied for beta.  Only the first value will be used.')
    beta <- beta[1]
  }

  betavals <- rep(beta, length(localities))
  names(betavals) <- paste0('beta_',localities)

  c(betavals, cpvals)
}

