library(CovMitigation)
library(ggplot2)
library(dplyr)

modeldir <- 'analysis/filter-updates.2020-10-25'
jan01 <- as.Date('2020-01-01')
months <- seq(as.Date('2020-11-01'), as.Date('2021-10-01'), by='month')
times <- as.numeric(months - jan01)
enddate <- as.numeric(as.Date('2021-01-31') - jan01)

## Scenario adjustments from the VA-Rand report
## NB: I've moved the vaccine timeline back by 1 month relative to what Rand is
##     suggesting.  I feel like December tends to be a lost month for many companies,
##     so I'd be surprised if a vaccine were rolled out starting in January.
##                N      D     J     F    M      A     M     J     J     A    S      O
seasonaladj <- c(1.10, 1.15, 1.15, 1.15, 1.10, 1.05, 1.00, 0.95, 0.85, 0.85, 0.95, 1.05)
traveladj <-   c(1.10, 1.15, 1.05, 1.00, 1.00, 1.00, 1.00, 1.00, 1.05, 1.05, 1.00, 1.00)
vaccineadj <-  c(1.00, 1.00, 1.00, 0.95, 0.95, 0.95, 0.90, 0.90, 0.90, 0.85, 0.85, 0.85)

bau_scen <- NULL
seasonality_scen <- Scenario(data.frame(locality=NA_character_, 
                                        time=times,
                                        parm='beta', isrelative=TRUE,
                                        value=seasonaladj))

travel_scen <- Scenario(data.frame(locality=NA_character_, 
                                        time=times,
                                        parm='beta', isrelative=TRUE,
                                        value=traveladj))

season_trav_scen <- Scenario(data.frame(locality=NA_character_, 
                                        time=times,
                                        parm='beta', isrelative=TRUE,
                                        value=seasonaladj*traveladj))

surge_scen <- Scenario(data.frame(locality=NA_character_,
                                  time=times[1],
                                  parm='beta', isrelative=TRUE,
                                  value=1.2))


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
