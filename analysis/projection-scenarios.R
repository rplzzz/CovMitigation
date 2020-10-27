library(CovMitigation)
library(ggplot2)
library(dplyr)

modeldir <- 'analysis/filter-updates.2020-10-18'
scenariodate <- 288
enddate <- 366

bau_scen <- NULL
statewide_surge_scen <- Scenario(data.frame(locality=NA_character_, time=scenariodate,
                                            parm='beta', isrelative=TRUE,
                                            value=1.2))
uva_surge_scen <-
    Scenario(dplyr::bind_rows(
        data.frame(locality='Charlottesvillecity', time=scenariodate, parm='beta',
                   isrelative=TRUE, value=1.2, stringsAsFactors=FALSE),
        data.frame(locality='AlbemarleCounty', time=scenariodate, parm='beta',
                   isrelative=TRUE, value=1.15, stringsAsFactors=FALSE)))


bau <- census_model_output(modeldir, enddate, 'BAU', scenario=bau_scen)
statewide_surge <- census_model_output(modeldir, enddate,
                                       'Statewide Surge 20%',
                                       scenario=statewide_surge_scen)
uva_surge <- census_model_output(modeldir, enddate, 'CHO-Alb Surge',
                                 scenario=uva_surge_scen)

combin_scenario_rslts <- dplyr::bind_rows(bau, statewide_surge, uva_surge)
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
