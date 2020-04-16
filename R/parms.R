
## Counties that seem to have a higher growth rate than the rest of the state
high_growth_counties <- 
  c('FairfaxCounty', 'ArlingtonCounty', 'PrinceWilliamCounty', 'HenricoCounty')

#' Model parameters and hyperparameters
#' 
#' @format 
#' \code{default_hparms} : Named list of default values for prior and likelihood
#' hyperparameters
#' @name parameters
#' @export
default_hparms <- list(fsalpha=100, fsbeta=100,
                       t0mulog=4.5, t0siglog=2,
                       bmulog=4, bsiglog=0.5,
                       nhosp_weight = 100,
                       nhosp_alpha = 12,
                       nhosp_beta = 190,
                       hg_counties = high_growth_counties)

#' Fill a vector with default values.
#' 
#' Add default values for any values not given in the first argument.
#' Works on vectors and lists.  Returns a structure of the same type as its first input.
#' 
#' @param parms Named vector of parameters
#' @param defaults Named vector of default parameters
#' @rdname parameters
#' @export
fill_defaults <- function(parms, defaults)
{
  for(n in names(defaults)) {
    if(! n %in% names(parms)) {
      parms[[n]] <- defaults[[n]]
    }
  }
  parms
}