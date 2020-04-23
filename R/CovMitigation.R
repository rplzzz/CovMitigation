
#' CovMitigation:  Data and modeling for UVAHS epidemic response
#' 
#' @section Modeling:
#' 
#' The main modeling function is \code{\link{run_scenario}}.  This function takes a
#' vector of time values(in days) to output and a list of parameters.
#' It produces a dataset giving key results,
#' including new infections and acute and ICU admissions.  Each call to the function runs a
#' single scenario, which is tagged with a name in the \code{scenario} column.  Therefore, 
#' concatenating the data frames from several runs and plotting with the color aesthetic 
#' set to "scenario" is an easy way to compare several different parameter runs.
#' 
#' The second important function is \code{\link{gen_post}}, which generates a function
#' that takes a vector of model parameters, runs the model and compares to observed data,
#' computes a log-likelihood, and adds a log-prior to get the log-posterior probability
#' density.  Embedded in this function is a secondary model for generating the expected
#' number of confirmed cases, given the model outputs.  Among other things, this 
#' calculation requires us to estimate the effect of targeted testing in raising the
#' number of expected cases, relative to what would be expected from a random sample.
#' 
#' 
#' @section Data:
#' 
#' Most of the package data is preprocessed data used in the modeling.  This includes
#' information such as demographic data for the local jurisdictions in Virginia,
#' historical UVAHS market shares by local jurisdiction, and so forth.
#' 
#' COVID-19 specific data is hosted in external packages.  We use the county-level
#' incidence reports from the New York Times, provided by the \code{NYTimesCOVID19}
#' package, and we use the state level testing counts from the COVID Tracking Project,
#' provided by the \code{vacovdata} package.  Further information about those datasets
#' can be found in the package documentation for those packages.
#' 
#' @docType package
#' @name CovMitigation
NULL
