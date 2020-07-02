library(CovMitigation)
library(tidyr)
library(dplyr)
library(ggplot2)
library(metrosamp)
library(doParallel)

registerDoParallel(8)

viscounties <- c('AlbemarleCounty', 'Charlottesvillecity', 'NelsonCounty',  ## Thomas Jefferson district
                 'FairfaxCounty', 'ArlingtonCounty', 'PrinceWilliamCounty',  ## Northern Virginia
                 'Richmondcity', 'HenricoCounty', 'ChesterfieldCounty'      ## Richmond area
                 )
  
day0fix <- 30

p0 <- c(eta=-0.7, xi=0, D0=7, A0=5, day_zero=day0fix, b=20, I0=10, Ts=4)

lpost <- gen_post(fixed_parms=p0)   # fixed parameters can still be overriden by supplying them explicitly

print(lpost(p0))

pstrt <- p0[-c(5)]
ctrl <- list(fnscale=-1, maxit=2000)
opt_uncons <- optim(pstrt, lpost, control=ctrl)
puncons <- opt_uncons$par

pmat <- rbind(puncons)
scen_names <- c('Unconstrained')

proj <- plt_projections(pmat, scen_names, p0, usedate=TRUE) + scale_color_brewer(type='qual')
print(proj)

modobs <- plt_modobs(pmat, scen_names, viscounties) + ggplot2::scale_color_brewer(type='qual')
print(modobs)
