## Diagnostic plots for model runs


#' Plot model predictions against data observations
#'
#' Plot model outputs and observed data.
#'
#' @param parms Parameters to run model for
#' @param scenarios Names of scenarios
#' @param counties Counties to include in the plot; if not specified, the top 4
#' by number of infections will be plotted.
#' @param maxdate Maximum date to include in plots.  Can be either a character
#'   string or a date object.
#' @param cols Number of columns in the grid of miniplots.
#' @param default_parms Default values to use for parameters not specified in parms
#' @importFrom dplyr %>%
#' @export
plt_modobs <- function(parms, scenarios=NULL, counties=NULL,
                       maxdate = NULL,
                       cols=3, default_parms=NULL)
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
  if(!is.null(maxdate)) {
    if(!inherits(maxdate, 'Date')) {
      maxdate <- as.Date(maxdate)
    }
    obsdata <- obsdata[obsdata$date <= maxdate, ]
  }

  wkdates <- unique(obsdata$date)

  modrslts <-
    lapply(seq_along(scenarios),
      function(i) {
        pvals <- fill_defaults(parms[i,], obs[['default_parm_vals']])
        modpvals <- as.list(pvals[! names(pvals) %in% c('b0', 'b1')])

        ## get output for every day up to the last in the dataset.
        tmax <- max(obsdata[['time']])
        tvals <- seq(infection_t0, tmax)
        modout <- run_scenario(tvals, modpvals, counties=counties)
        modout[['scenario']] <- scenarios[i]
        modout[['fi']] <- modout[['PopInfection']] / modout[['population']]
        modout[['date']] <- modout[['time']] + as.Date('2020-01-01')

        modout[['b0']] <- pvals[['b0']]
        modout[['b1']] <- pvals[['b1']]

        modout <- modout[,c('scenario','date','time','locality','fi', 'b0', 'b1')]

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

  pltdata[['bias']] <- pltdata[['b0']] - pltdata[['b1']] * log(pltdata[['ntest']])

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
    newSympto <- population <- popTotal <- scenario <- value <- NULL

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
               modparms <- as.list(pset[! names(pset) %in% c('b0', 'b1')])
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
                       popTotal = sum(population*marketFraction)) %>%
      dplyr::mutate(PopPrevalence = PopInfection / popTotal)
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
                       popTotal = sum(population)) %>%
      dplyr::mutate(PopPrevalence = PopInfection / popTotal)
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


#' Generate diagnostic plots for filter models
#'
#' @param modelfit Modelfit element of the filter-fit structure
#' @param obsrun Optional data frame or matrix of results from the run that generated
#' the data.
#' @return A list of ggplot objects
#' @export
filter_model_diagnositc_plots <- function(modelfit, obsrun=NULL)
{
  prevplot <-
    ggplot2::ggplot(modelfit, ggplot2::aes(x=time)) +
    ggplot2::geom_line(ggplot2::aes(y=fi, color='model fit'), size=1.2) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin=filo, ymax=fihi), alpha=0.5) +
    ggplot2::ylab('prevalence') +
    ggplot2::theme_bw()

  betaplot <-
    ggplot2::ggplot(modelfit, ggplot2::aes(x=time)) +
    ggplot2::geom_line(ggplot2::aes(y=beta), size=1.2) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin=betalo, ymax=betahi), alpha=0.5) +
    ggplot2::ylab('beta') +
    ggplot2::theme_bw()

  importplot <-
    ggplot2::ggplot(modelfit, ggplot2::aes(x=time)) +
    ggplot2::geom_line(ggplot2::aes(y=import_cases), size=1.2) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin=import_caseslo, ymax=import_caseshi),
                         alpha=0.5) +
    ggplot2::ylab('Import cases') +
    ggplot2::theme_bw()

  if(!is.null(obsrun)) {
    if(!is.data.frame(obsrun)) {
      stopifnot(is.matrix(obsrun))
      obsrun <- as.data.frame(obsrun)
    }
    totpop <- obsrun$S + obsrun$E + obsrun$I + obsrun$Is + obsrun$R
    obsrun$fi <- (obsrun$I + obsrun$Is) / totpop
    prevplot <- prevplot +
      ggplot2::geom_line(data=obsrun, ggplot2::aes(y=fi, color='simulated run'), size=1.2)
  }

  list(prevalence=prevplot, beta=betaplot, import=importplot)
}


#' Plot uncertainty ranges by vintage
#'
#' @param vintagedata Output from \code{\link{project_filter_model}}.
#' @param what Variable to plot.  S, I, Is, Itot, and fi are available.
#' @export
vintage_plot <- function(vintagedata, what)
{
  nvintmax <- 5

  locol <- paste0(what, 'lo')
  hicol <- paste0(what, 'hi')

  pltdata <- vintagedata[ , c('vintage','time', what, locol, hicol)]
  names(pltdata) <- c('vintage','time','y', 'ylo', 'yhi')

  vintages <- unique(pltdata[['vintage']])
  ## If there are too many to plot easily, pick a subset as evenly distributed as possible.
  if(length(vintages) > nvintmax) {
    vintages <- vintages[round(seq(1, length(vintages), length.out=nvintmax))]
    pltdata <- pltdata[pltdata[['vintage']] %in% vintages, ]
  }
  pltdata[['vintage']] <- as.factor(pltdata[['vintage']])

  ggplot2::ggplot(pltdata, ggplot2::aes(x=time, color=vintage, fill=vintage)) +
    ggplot2::geom_line(mapping=ggplot2::aes(y=y), size=1.2) +
    ggplot2::geom_ribbon(mapping=ggplot2::aes(ymin=ylo, ymax=yhi), alpha=0.25) +
    ggplot2::theme_bw() +
    ggplot2::ylab(what) +
    ggthemes::scale_color_solarized() +
    ggthemes::scale_fill_solarized()
}
