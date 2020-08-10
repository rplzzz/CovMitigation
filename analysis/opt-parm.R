library(CovMitigation)
library(tidyr)
library(dplyr)
library(ggplot2)
library(metrosamp)
library(doParallel)

registerDoParallel(6)

viscounties <- c('AlbemarleCounty', 'Charlottesvillecity', 'NelsonCounty',  ## Thomas Jefferson district
                 'FairfaxCounty', 'ArlingtonCounty', 'PrinceWilliamCounty',  ## Northern Virginia
                 'Richmondcity', 'HenricoCounty', 'ChesterfieldCounty',      ## Richmond area
                 'VirginiaBeachcity', 'Portsmouthcity', 'Norfolkcity'       ## SE VA
)

cpvals <- c(D0=5.5, A0=2, Ts=5.5, mask_effect=-0.5, b0=100, b1=1)
p0 <- gen_parm_vec(common_parms = cpvals, beta = 0.2)
## Raise beta in some of the statewide hotspots
chgparms <- paste0('beta_', 
                   c('FairfaxCounty', 'Fairfaxcity','FallsChurchcity', 'ArlingtonCounty', 
                     'Alexandriacity',
                     'PrinceWilliamCounty', 'Manassascity', 'ManassasParkcity',
                     'Richmondcity', 'HenricoCounty', 'ChesterfieldCounty'))
p0[chgparms] <- 0.25

lpost <- gen_post(fixed_parms=p0, maxdate = as.Date('2020-06-30'))

print(lpost(p0))

pstrt <- p0
ctrl <- list(fnscale=-1, maxit=2000)
opt_uncons <- optim(pstrt, lpost, control=ctrl)
puncons <- opt_uncons$par

opt_uncons2 <- optim(puncons, lpost, control=list(fnscale=-1, maxit=500), method='BFGS')
puncons2 <- opt_uncons$par

