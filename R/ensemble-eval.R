#### Functions for evaluating ensembles of covid models

### Helper function for ensemble eval:  get observations and predictions by time
### for each ensemble, for a single locality.
get_ldata <- function(m) 
{
  ensemble_fcst <- m$history[ , c('time', 'id', 'ncase')]
  obsdata <- m$obsdata[ , c('time', 'fips', 'ntest', 'npos')]
  obsdata$pop <- vdhcovid::getpop(fips=obsdata$fips)
  rslt <- dplyr::left_join(ensemble_fcst, obsdata, by='time')
  if(!all(complete.cases(rslt))) {
    warning('get_ldata:  Some times missing observations for locality: ', m$obsdata$locality[1])
  }
  rslt[rslt$ntest > 0, ]
}


#' Evaluate ensemble members according to their agreement with observations.
#' 
#' For each ensemble member, compare its predicted number of cases with the
#' observed number of cases and compute a log-likelihood function.
#' 
#' Right now we're just looking at number of cases, given the effective number
#' of tests.  Eventually we hope to include the number of UVA hospitalizations
#' as a constraint.
#' 
#' @param modlist List of \code{filter-fit} objects. Each fit is for a single 
#' locality.
#' @return Table of ensemble ID numbers and log-likelihood values, along with some
#' diagnostic statistics.
#' @export
ensemble_eval <- function(modlist)
{
  fips <- logl <- id <- meanlogl <- varlogl <- NULL  # silence warnings
  
  cmpdata <- dplyr::bind_rows(lapply(modlist, get_ldata))
  p <- cmpdata$ncase / cmpdata$ntest
  m <- ceiling(p*cmpdata$pop)
  n <- cmpdata$pop-m
  k <- cmpdata$ntest
  x <- cmpdata$npos
  cmpdata$logl <- dhyper(x, m, n, k, log=TRUE)
  
  ## Aggregate to the county level first so that we can produce some statistics
  ## on whether or not the aggregate log likelihoods are being driven by outlier
  ## counties.
  intermeddata <- 
    dplyr::group_by(cmpdata, fips, id) %>%
    dplyr::summarise(logl = sum(logl)) %>%
    dplyr::ungroup()
  
  dplyr::group_by(intermeddata, id) %>%
    dplyr::summarise(meanlogl = mean(logl),
                     varlogl = var(logl),
                     zm2 = (as.numeric(quantile(logl, pnorm(-2))) - meanlogl) / varlogl,
                     zp2 = (as.numeric(quantile(logl, pnorm(2))) - meanlogl) / varlogl,
                     logl = sum(logl)
                     )
}
