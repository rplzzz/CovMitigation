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
  idx <- match(localities, vdhcovid::valocalities$locality)
  
  exp(parms[['eta']] + parms[['xi']]*vdhcovid::valocalities[['stdpopdens']][idx])
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
