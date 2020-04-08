library(magrittr)
library(CovMitigation)
library(doParallel)

T0vals <- seq(4.5,7.5)
registerDoParallel(8)
baseparms <- list()

basescens <- 
  dplyr::bind_rows(
    lapply(T0vals, function(T0) {
      parms <- c(baseparms, list(T0=T0))
      run_scenario(270, parms, paste('td = ', T0))
    }))

pltdata <- basescens %>%
  dplyr::group_by(time, scenario) %>%
  dplyr::summarise(
    newCases = sum(newCases*marketFraction),
    PopSympto = sum(PopSympto*marketFraction),
    PopInfection = sum(PopInfection*marketFraction),
    PopCumulInfection = sum(PopCumulInfection*marketFraction),
    PopCumulFrac = sum(PopCumulInfection)/sum(population*marketFraction),
    popTotal = sum(population*marketFraction)
  )
