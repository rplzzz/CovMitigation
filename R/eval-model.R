#' Compute the deviance for models in an ensemble.
#'
#' The deviance for a model is 2(log L - log L0), where L is the likelihood
#' value for the model and L0 is the likelihood for the saturated model (a
#' hypothetical model that gets every observation exactly correct).
#'
#' Theoretically, we should compute the expected observations by averaging the
#' prevalence over the week corresponding to the observation, but we get a
#' pretty good approximation to that by just taking the point estimates at the
#' end of the week.
#'
#' @param filterfit Filter-fit object for the ensemble.  The history must be
#'   included in order to calculate the deviance (see
#'   \code{\link{fit_filter}}).
#' @return The input object, with a new entry \code{deviance} added.
#' @export
ensemble_deviance <- function(filterfit)
{
  obsdata <- filterfit[['obsdata']]
  history <- filterfit[['history']]

  history <- history[history[['time']] %in% obsdata[['time']], ]

  modata <- dplyr::left_join(history, obsdata[ , c('time','ntest','npos')], by='time')
  ids <- filterfit[['ids']]
  fs <- filterfit[['finalstate']]
  b0 <- fs[ ,'b0']
  b1 <- fs[ ,'b1']

  modata <- dplyr::left_join(modata, data.frame(id=ids, b0=b0, b1=b1), by='id')
  b <- pmax(1, modata[['b0']] - modata[['b1']] * log(modata[['ntest']]))

  fibiased <- padjust(modata[['fi']], b)
  fpos <- modata[['npos']] / modata[['ntest']]

  d0 <- dbinom(modata[['npos']], modata[['ntest']], fpos, log=TRUE)
  modata[['dev']] <- 2*(dbinom(modata[['npos']], modata[['ntest']], fibiased,
                               log=TRUE) - d0)

  filterfit[['deviance']] <-
      dplyr::group_by(modata, id) %>%
      dplyr::summarise(deviance=sum(dev)) %>%
      dplyr::arrange(dplyr::desc(deviance))

  filterfit
}
