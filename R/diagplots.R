## Diagnostic plots for model runs

fill_defaults <- function(parms, defaults)
{
  
}

#' Plot model predictions against data observations
#' 
#' Plot model outputs and observed data.
#' 
#' @param parms Parameters to run model for
#' @param obs Observed data as returned by \code{\link{get_obsdata}}
#' @param ncounty Number of counties to include in the plot
#' @importFrom dplyr %>%
#' @export
plt_modobs <- function(parms, ncounty=4)
{
  ## silence warnings
  fips <- newcases <- predicted <- mtype <- otype <- NULL
  
  obs <- get_obsdata()
  obsdata <- obs[['obsdata']]
  modparms <- as.list(parms[! names(parms) %in% c('day_zero', 'b')])
  
  if('day_zero' %in% names(parms)) {
    day0 <- parms[['day_zero']]
  }  
  else {
    day0 <- obs[['default_parm_vals']][['day_zero']]
  }
  obsdata[['time']] <- obsdata[['day']] - day0
  
  if('b' %in% names(parms)) {
    b <- parms[['b']]
  }
  else {
    b <- obs[['default_parm_vals']][['b']]
  }
  
  ## get output for every day up to the last in the dataset.
  tmax <- max(obsdata[['time']])
  tvals <- c(0, seq(to = tmax, length.out = floor(tmax)))
  
  modout <- run_scenario(tvals, modparms)
  
  obs[['obsdata']] <- obsdata
  
  cmp <- align_modout(modout, obs)

  ncases <- dplyr::group_by(cmp, fips) %>%
    dplyr::summarise(ncases = sum(newcases)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(dplyr::desc(ncases))
  
  use_fips <- ncases[['fips']][seq(1,ncounty)]
  
  pltdata <- cmp[cmp[['fips']] %in% use_fips,]
  
  pltdata[['predicted']] <- pltdata[['model.newcases']] * pltdata[['ftest']] * b
  pltdata[['otype']] <- "observed"
  pltdata[['mtype']] <- "predicted"
  
  ggplot2::ggplot(data=pltdata, ggplot2::aes(x=date)) + 
    ggplot2::geom_line(mapping=ggplot2::aes(y=predicted, linetype=mtype), size=1.2) + 
    ggplot2::geom_point(mapping=ggplot2::aes(y=newcases, shape=otype)) + 
    ggplot2::facet_wrap(~Locality) +
    ggplot2::guides(shape=ggplot2::guide_legend(''), linetype=ggplot2::guide_legend('')) +
    ggplot2::ylab('New Cases') +
    ggplot2::theme_bw()
}

#' Plot model projections for comparison
#' 
#' All of the results we plot here are adjusted for the UVA market fraction.
#' 
#' @param parms Matrix of combinations of parameters.
#' @param scenarios Vector of scenario names
#' @param tmax Maximum time to run to
#' @param what What output to plot.  Options are 
#' newCases, newSympto, PopSympto, PopInfection, PopCumulInfection, PopCumulFrac, popTotal
#' @export
plt_projections <- function(parms, scenarios, tmax=270, what='newSympto')
{
  ## silence warnings
  newCases <- marketFraction <- PopSympto <- PopInfection <- PopCumulInfection <-
    population <- scenario <- value <- NULL
  
  if(!is.matrix(parms)) {
    parms <- matrix(parms, nrow=1)
  }
  
  stopifnot(nrow(parms) == length(scenarios))
  
  tvals <- seq(0, tmax)
  
  modouts <- 
    dplyr::bind_rows(
      lapply(seq(1, nrow(parms)), 
             function(i) {
               pset <- parms[i,]
               modparms <- as.list(pset[! names(pset) %in% c('day_zero', 'b')])
               modout <- run_scenario(tvals, modparms, scenarios[i])
               if(!is.na(pset['day_zero'])) {
                 modout[['time']] <- modout[['time']] + pset['day_zero']
               }
               modout
             })  
    )

  modouts <- dplyr::group_by(modouts, time, scenario)
  pltdata <-
    dplyr::summarise(modouts,
                     newCases = sum(newCases*marketFraction),
                     newSympto = sum(newSympto*marketFraction),
                     PopSympto = sum(PopSympto*marketFraction),
                     PopInfection = sum(PopInfection*marketFraction),
                     PopCumulInfection = sum(PopCumulInfection*marketFraction),
                     PopCumulFrac = sum(PopCumulInfection)/sum(population*marketFraction),
                     popTotal = sum(population*marketFraction)
    )
  
  pltdata[['value']] <- pltdata[[what]]
  
  ggplot2::ggplot(data=pltdata, ggplot2::aes(x=time, y=value, color=scenario)) +
    ggplot2::geom_line(size=1.2) +
    ggplot2::xlab(what) +
    ggplot2::theme_bw()
}