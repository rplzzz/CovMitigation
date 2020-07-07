#' Indicator function for mask use requirement
#' 
#' For an input vector of dates (or times given as days since 2020-01-01), return 
#' 1 for dates on which mask use was required, zero for dates on which it was not.
#' 
#' Executive order 63 established a requirement to wear masks in most public 
#' places in the Commonwealth of Virginia, starting Friday, 29 May 2020.
#' 
#' @param t Vector of dates or times (days since Jan 01)
#' @export
mask_indicator <- function(t) 
{
  UseMethod('mask_indicator', t)
}

#' @describeIn mask_indicator Default version
#' @export
mask_indicator.default <- function(t) 
{
  as.numeric(t >= 149)       # =1 for dates on or after 29 May 2020.
}

#' @describeIn mask_indicator Date class version
#' @export
mask_indicator.Date <- function(t) 
{
  as.numeric(t >= as.Date('2020-05-29'))
}
