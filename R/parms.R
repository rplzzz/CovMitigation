
#' Model parameters and hyperparameters
#' 
#' @format 
#' \code{default_hparms} : Named list of default values for prior and likelihood
#' hyperparameters
#' @name parameters
#' @export
default_hparms <- list(fsalpha=100, fsbeta=100,
                       t0mulog=log(c(5,10,15,40)), t0siglog=c(1, 1.5, 1.5, 2),
                       bmulog=4, bsiglog=0.5,
                       nhosp_weight = 100,
                       nhosp_alpha = 12,
                       nhosp_beta = 190)

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
