library(CovMitigation)
library(ggplot2)

bau_scen <- NULL
statewide_surge_scen <- Scenario(data.frame(locality=NA_character_, time=251,
                                            parm='beta', isrelative=TRUE,
                                            value=1.25))
uva_surge_scen <-
    Scenario(dplyr::bind_rows(
        data.frame(locality='Charlottesvillecity', time=250, parm='beta',
                   isrelative=TRUE, value=1.25, stringsAsFactors=FALSE),
        data.frame(locality='Charlottesvillecity', time=250, parm='gamma',
                   isrelative=TRUE, value=1/0.9, stringsAsFactors = FALSE),
        data.frame(locality='AlbemarleCounty', time=250, parm='beta',
                   isrelative=TRUE, value=1.15, stringsAsFactors=FALSE)))


bau <- census_model_output('analysis/filter-fits', 305, 'BAU', scenario=bau_scen)
statewide_surge <- census_model_output('analysis/filter-fits', 305,
                                       'Statewide Surge 25%',
                                       scenario=statewide_surge_scen)
uva_surge <- census_model_output('analysis/filter-fits', 305, 'CHO-Alb Surge',
                                 scenario=uva_surge_scen)

combin_scenario_rslts <- dplyr::bind_rows(bau, statewide_surge, uva_surge)
combin_plt <- dplyr::group_by(combin_scenario_rslts, time, date, scenario) %>%
  dplyr::summarise(PopInfection = median(PopInfection), PopSympto = median(PopSympto),
                   fracInfection = median(fracInfection), newCases = median(newCases),
                   newSympto = median(newSympto))

plt_ns <- ggplot(combin_plt, aes(x=date, y=newSympto, color=scenario)) +
  geom_line(size=1.2) +
  ggtitle('UVA Market Share Adjusted') +
  ylab('Daily New Symptomatic Cases') +
  theme_bw()

plt_prev <- ggplot(combin_plt, aes(x=date, y=fracInfection, color=scenario)) +
  geom_line(size=1.2) +
  ggtitle('UVA Market Share Adjusted') +
  ylab('Prevalence') +
  theme_bw()

print(plt_ns)
print(plt_prev)
