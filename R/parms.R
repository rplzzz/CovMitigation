#' Model parameters and hyperparameters
#' 
#' @format 
#' \code{default_hparms} : Named list of default values for prior and likelihood
#' hyperparameters
#' @name parameters
#' @export
default_hparms <- list(fsalpha=100, fsbeta=100,
                       t0shp=14, t0rate=2,
                       bmulog=4, bsiglog=0.5,
                       nhosp_weight = 100,
                       nhosp_alpha = 4,
                       nhosp_beta = 98)

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
