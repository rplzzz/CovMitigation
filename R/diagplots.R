## Diagnostic plots for model runs


#' Plot model predictions against data observations
#' 
#' Plot model outputs and observed data.
#' 
#' @param parms Parameters to run model for
#' @param scenarios Names of scenarios
#' @param counties Counties to include in the plot; if not specified, the top 4
#' by number of infections will be plotted.
#' @param cols Number of columns in the grid of miniplots.
#' @param default_parms Default values to use for parameters not specified in parms
#' @importFrom dplyr %>%
#' @export
plt_modobs <- function(parms, scenarios=NULL, counties=NULL, cols=3, default_parms=NULL)
{
  ## silence warnings
  fips <- newcases <- predicted <- NULL
  
  obs <- get_obsdata()
  if(!is.null(default_parms)) {
    obs$default_parm_vals <- fill_defaults(default_parms, obs$default_parms_vals)
  }
  if(!is.matrix(parms)) {
    parms <- t(as.matrix(parms))
  }
  if(is.null(scenarios)) {
    scenarios <- paste0('s', seq(1,nrow(parms)))
  }

    
  
  obsdata <- obs[['obsdata']]
  wkdates <- unique(obsdata$date)
  
  modrslts <- 
    lapply(seq_along(scenarios),
      function(i) {
        pvals <- fill_defaults(parms[i,], obs[['default_parm_vals']])
        modpvals <- as.list(pvals[! names(pvals) %in% c('b')])
        
        b <- pvals[['b']]
        ## get output for every day up to the last in the dataset.
        tmax <- max(obsdata[['time']])
        tvals <- seq(infection_t0, tmax)
        modout <- run_scenario(tvals, modpvals, counties=counties)
        modout[['scenario']] <- scenarios[i]
        modout[['bias']] <- pvals[['b']]
        modout[['fi']] <- modout[['PopInfection']] / modout[['population']]
        modout[['date']] <- modout[['time']] + as.Date('2020-01-01')
      
        modout <- modout[,c('scenario','date','time','locality','fi','bias')]
        
        lmout <- split(modout, modout$locality, drop=TRUE)
        
        dplyr::bind_rows(
          lapply(lmout, function(df) {wkagg(wkdates, df, 'fi')})
        )
        
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
  
  pltdata <- 
    dplyr::bind_rows(modrslts) %>%
    dplyr::left_join(obsdata, by=c('date','time', 'locality'))
  pltdata <- pltdata[pltdata[['fips']] %in% use_fips,]
  pltdata <- pltdata[!is.na(pltdata[['npos']]), ]
  
  pltdata[['predicted']] <- 
    padjust(pltdata[['fi']], pltdata[['bias']]) * pltdata[['ntest']]
  
  pltdata[['locality']] <- ordered(pltdata[['locality']], levels = counties)
  
  ggplot2::ggplot(data=pltdata, ggplot2::aes(x=date)) + 
    ggplot2::geom_line(mapping=ggplot2::aes(y=predicted, linetype='predicted', color=scenario), size=1.2) + 
    ggplot2::geom_point(mapping=ggplot2::aes(y=npos, shape='observed')) + 
    ggplot2::facet_wrap(~locality, scales='free_y', ncol=cols) +
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
#' @param default_parms Default values for parameters not mentioned in parms
#' @param tmax Maximum time to run to
#' @param what What output to plot.  Options are 
#' newCases, newSympto, PopSympto, PopInfection, PopCumulInfection, PopCumulFrac, popTotal
#' @param usedate If \code{TRUE}, display date on the x-axis; otherwise, display model
#' time.
#' @param counties If specified, filter to just the requested counties.  Otherwise,
#' include the whole state.
#' @param marketadjust If \code{TRUE}, adjust projections for UVAHS market share; 
#' otherwise, show raw totals.
#' @export
plt_projections <- function(parms, scenarios, default_parms = list(), tmax=270, what='newSympto', usedate=TRUE, 
                            counties=NULL, marketadjust=FALSE)
{
  ## silence warnings
  newCases <- marketFraction <- PopSympto <- PopInfection <- PopCumulInfection <-
    newSympto <- population <- scenario <- value <- NULL
  
  if(!is.matrix(parms)) {
    parms <- t(as.matrix(parms))
  }
  
  stopifnot(nrow(parms) == length(scenarios))
  
  tvals <- seq(0, tmax)
  
  modouts <- 
    dplyr::bind_rows(
      lapply(seq(1, nrow(parms)), 
             function(i) {
               pset <- fill_defaults(parms[i,], default_parms)
               modparms <- as.list(pset[! names(pset) %in% c('b')])
               modout <- run_scenario(tvals, modparms, scenarioName = scenarios[i])
               if(!is.null(counties)) {
                 modout <- modout[modout$locality %in% counties,]
               }
               modout
             })  
    )


  modouts <- dplyr::group_by(modouts, time, scenario)
    
  if(marketadjust) {
    pltdata <-
      dplyr::summarise(modouts,
                       newCases = sum(newCases*marketFraction),
                       newSympto = sum(newSympto*marketFraction),
                       PopSympto = sum(PopSympto*marketFraction),
                       PopInfection = sum(PopInfection*marketFraction),
                       PopCumulInfection = sum(PopCumulInfection*marketFraction),
                       PopCumulFrac = sum(PopCumulInfection)/sum(population*marketFraction),
                       popTotal = sum(population*marketFraction))
  }
  else {
    pltdata <-
      dplyr::summarise(modouts,
                       newCases = sum(newCases),
                       newSympto = sum(newSympto),
                       PopSympto = sum(PopSympto),
                       PopInfection = sum(PopInfection),
                       PopCumulInfection = sum(PopCumulInfection),
                       PopCumulFrac = sum(PopCumulInfection)/sum(population),
                       popTotal = sum(population))
    
  }
  
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
