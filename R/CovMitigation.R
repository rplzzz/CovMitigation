
#' CovMitigation:  Data and modeling for UVAHS epidemic response
#' 
#' @section Modeling:
#' 
#' The main modeling function is \code{\link{run_scenario}}.  This function takes an
#' end time (in days) and a list of parameters and produces a dataset giving key results,
#' including new infections and acute and ICU admissions.  Each call to the function runs a
#' single scenario, which is tagged with a name in the \code{scenario} column.  Therefore, 
#' concatenating the data frames from several runs and plotting with the color aesthetic 
#' set to "scenario" is an easy way to compare several different parameter runs.
#' 
#' @section Data:
#' 
#' Most of the package data is preprocessed data used in the modeling; however, the 
#' \code{\link{vdh_covid19}} dataset contains daily data (starting from 25 March) from
#' the Virginia department of health on the outbreak.  Since this dataset is frequently
#' updated, we don't necessarily generate a release or increment the package version 
#' every time we update it.
#' 
#' @docType package
#' @name CovMitigation
NULL
