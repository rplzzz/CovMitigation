## Diagnostic plots for model runs


#' Plot model predictions against data observations
#' 
#' Plot model outputs and observed data.
#' 
#' @param parms Parameters to run model for
#' @param scenarios Names of scenarios
#' @param counties Counties to include in the plot; if not specified, the top 4
#' by number of infections will be plotted.
#' @param default_parms Default values to use for parameters not specified in parms
#' @importFrom dplyr %>%
#' @export
plt_modobs <- function(parms, scenarios=NULL, counties=NULL, default_parms=NULL)
{
  ## silence warnings
  fips <- newcases <- predicted <- NULL
  
  obs <- get_obsdata()
  if(!is.null(default_parms)) {
    obs$default_parm_vals <- fill_defaults(obs$default_parm_vals, default_parms)
  }
  if(!is.matrix(parms)) {
    parms <- t(as.matrix(parms))
  }
  if(is.null(scenarios)) {
    scenarios <- paste0('s', seq(1,nrow(parms)))
  }
  
  
  obsdata <- obs[['obsdata']]
  
  modrslts <- 
    lapply(seq_along(scenarios),
      function(i) {
        pvals <- fill_defaults(parms[i,], obs[['default_parm_vals']])
        modpvals <- as.list(pvals[! names(pvals) %in% c('day_zero', 'b')])
  
        day0 <- pvals[['day_zero']]
        b <- pvals[['b']]
        ## get output for every day up to the last in the dataset.
        tmax <- max(obsdata[['time']])
        tvals <- c(day0, seq(ceiling(day0), tmax))
        modout <- run_scenario(tvals, modpvals)
        modout[['scenario']] <- scenarios[i]
        modout[['bias']] <- pvals[['b']]
        modout
      })
  
  if (is.null(counties)){
    ncounty <- 4
    ncases <- dplyr::group_by(obsdata, fips) %>%
      dplyr::summarise(ncases = sum(newcases)) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(dplyr::desc(ncases))
    
    use_fips <- ncases[['fips']][seq(1,ncounty)]
  } 
  else {
    use_fips <- 
      va_county_first_case$FIPS[va_county_first_case$Locality %in% counties]
  }
  
  modrslts <- 
    dplyr::bind_rows(modrslts) %>%
    dplyr::left_join(obs[['fips_codes']], by=c(locality='Locality'))
  modrslts[['Itot']] <- modrslts[['I']] + modrslts[['Is']]
  modrslts[['fi']] <- modrslts[['Itot']] / modrslts[['population']]
  
  pltdata <- dplyr::left_join(modrslts, obsdata, by=c('fips','time'))
  pltdata <- pltdata[pltdata[['fips']] %in% use_fips,]
  pltdata <- pltdata[!is.na(pltdata[['newcases']]), ]
  pltdata[['popfrac']] <- pltdata[['population']] / vaMyAgeBands$Total[1]
  
  pltdata[['predicted']] <- 
    padjust(pltdata[['fi']], pltdata[['bias']]) * pltdata[['ntest']] * pltdata[['popfrac']]
  
  ggplot2::ggplot(data=pltdata, ggplot2::aes(x=date)) + 
    ggplot2::geom_line(mapping=ggplot2::aes(y=predicted, linetype='predicted', color=scenario), size=1.2) + 
    ggplot2::geom_point(mapping=ggplot2::aes(y=newcases, shape='observed')) + 
    ggplot2::facet_wrap(~locality) +
    #ggplot2::labs(color='Scenario', linetype='', shape='') +
    ggplot2::guides(
      color=ggplot2::guide_legend('Scenario', order=1),
      shape=ggplot2::guide_legend('',order=3), 
      linetype=ggplot2::guide_legend('',order=2)) +
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
#' @param usedate If TRUE, reference model time to 1-Jan
#' @export
plt_projections <- function(parms, scenarios, tmax=270, what='newSympto', usedate=FALSE)
{
  ## silence warnings
  newCases <- marketFraction <- PopSympto <- PopInfection <- PopCumulInfection <-
    population <- scenario <- value <- NULL
  
  if(!is.matrix(parms)) {
    parms <- t(as.matrix(parms))
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
  pltdata[['date']] <- pltdata[['time']] + as.Date('2020-01-01')
  
  if(usedate) {
    ggplot2::ggplot(data=pltdata, ggplot2::aes(x=date, y=value, color=scenario)) +
      ggplot2::geom_line(size=1.2) +
      ggplot2::ylab(what) +
      ggplot2::theme_bw()
  }
  else {
    ggplot2::ggplot(data=pltdata, ggplot2::aes(x=time, y=value, color=scenario)) +
      ggplot2::geom_line(size=1.2) +
      ggplot2::ylab(what) +
      ggplot2::theme_bw()
    }
}
