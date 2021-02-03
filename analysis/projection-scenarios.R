library(CovMitigation)
library(ggplot2)
library(dplyr)

modeldate <- '2021-01-24'
modeldir <- paste0('analysis/filter-updates.', modeldate)
message('model date: ', modeldate, '\tmodel dir: ', modeldir)
jan01 <- as.Date('2020-01-01')
nov01 <- as.Date('2020-11-01')
months <- seq(as.Date('2020-11-01'), as.Date('2021-10-01'), by='month')
times <- as.numeric(months - jan01)
enddate <- as.numeric(as.Date('2021-04-30') - jan01)
lastobsdate <- as.Date('2021-01-27')
modelmonth <- lubridate::as.period(lubridate::interval(nov01, as.Date(modeldate)))@month

## Start of the scenario adjustment
scenariodate <- as.numeric(lastobsdate - jan01) + 1
times <- pmax(times, scenariodate)

## Sigmoid function for vaccine adjustment
vadjust <- function(x1, x2, sat) {
  x0 <- 0.5*(x1 + x2)
  k <- log(100*sat - 1) / (x0-x1)
  function(x) {
    v <- ifelse(x >= x1, sat / (1 + exp(-k*(x-x0))), 0)
    1 - v
  }
}

## Scenario adjustments from the VA-Rand report
## NB: I've moved the vaccine timeline back by 1 month relative to what Rand is
##     suggesting.  I feel like December tends to be a lost month for many companies,
##     so I'd be surprised if a vaccine were rolled out starting in January.
##                N      D     J     F    M      A     M     J     J     A    S      O
seasonaladj <- c(1.10, 1.15, 1.15, 1.15, 1.10, 1.05, 1.00, 0.95, 0.85, 0.85, 0.95, 1.05)
traveladj <-   c(1.10, 1.15, 1.05, 1.00, 1.00, 1.00, 1.00, 1.00, 1.05, 1.05, 1.00, 1.00)
vaccineadj <- vadjust(4, 11, 0.75)(seq(1,12))

## Once the model has been running for a while, it should reflect these projections
## for seasonal, travel, and vaccine adjustments.  Therefore, we scale the adjustments
## according to what should be baked into the model already.
if(modelmonth >= 1 && modelmonth <= 12) {
  seasonaladj <- seasonaladj / seasonaladj[modelmonth]
  traveladj <- traveladj / traveladj[modelmonth]
  vaccineadj <- vaccineadj / vaccineadj[modelmonth]
}

## As of the 5 Nov. updates, it looks like the beta values in Charlottesville and
## Albemarle are being depressed by the ramp-up in UVA testing.  This adjustment
## offsets that effect
# bau_scen <-
#   Scenario(dplyr::bind_rows(
#     data.frame(locality='Charlottesvillecity', time=scenariodate, parm='beta',
#                isrelative=TRUE, value=1.2*vaccineadj, stringsAsFactors=FALSE),
#     data.frame(locality='AlbemarleCounty', time=scenariodate, parm='beta',
#                isrelative=TRUE, value=1.2*vaccineadj, stringsAsFactors=FALSE)))

bau_scen <- Scenario(data.frame(locality=NA_character_,
                     time=times,
                     parm='beta', isrelative=TRUE,
                     value=vaccineadj))

seasonality_scen <- Scenario(data.frame(locality=NA_character_, 
                                        time=times,
                                        parm='beta', isrelative=TRUE,
                                        value=seasonaladj*vaccineadj))

travel_scen <- Scenario(data.frame(locality=NA_character_, 
                                        time=times,
                                        parm='beta', isrelative=TRUE,
                                        value=traveladj*vaccineadj))

season_trav_scen <- Scenario(data.frame(locality=NA_character_, 
                                        time=times,
                                        parm='beta', isrelative=TRUE,
                                        value=seasonaladj*traveladj*vaccineadj))

surge_scen <- Scenario(data.frame(locality=NA_character_,
                                  time=times,
                                  parm='beta', isrelative=TRUE,
                                  value=1.2*vaccineadj))


bau <- census_model_output(modeldir, enddate, 'BAU', scenario=bau_scen)
season <- census_model_output(modeldir, enddate,
                              'Seasonality adjustment',
                              scenario=seasonality_scen)
travel <- census_model_output(modeldir, enddate, 'Travel adjustment',
                              scenario=travel_scen)
season_travel <- census_model_output(modeldir, enddate, 'Seasonality and Travel',
                                    scenario=season_trav_scen)
surge <- census_model_output(modeldir, enddate,
                             'Flat 20% Surge',
                             scenario=surge_scen)

combin_scenario_rslts <- dplyr::bind_rows(bau, season, travel, season_travel, surge)
combin_plt <- dplyr::group_by(combin_scenario_rslts, time, date, scenario) %>%
  dplyr::summarise(PopInfection = median(PopInfection), PopSympto = median(PopSympto),
                   fracInfection = median(fracInfection), newCases = median(newCases),
                   newSympto = median(newSympto))

newsympto_bounds <- c(0, max(combin_plt$newSympto))
prev_bounds <- c(0, max(combin_plt$fracInfection))

plt_ns <- ggplot(combin_plt, aes(x=date, y=newSympto, color=scenario)) +
  geom_line(size=1.2) +
  ggtitle('UVA Market Share Adjusted') +
  ylab('Daily New Symptomatic Cases') +
  coord_cartesian(ylim=newsympto_bounds) +
  theme_bw()

plt_prev <- ggplot(combin_plt, aes(x=date, y=fracInfection, color=scenario)) +
  geom_line(size=1.2) +
  ggtitle('UVA Market Share Adjusted') +
  ylab('Prevalence') +
  coord_cartesian(ylim=prev_bounds) +
  theme_bw()

print(plt_ns)
print(plt_prev)
