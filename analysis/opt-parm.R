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

p0 <- c(eta=-1.25, xi=0.1, zeta=0, D0=7, A0=5, I0=10, Ts=4, mask_effect=0,
        b0=20, b1=0.1)

lpost <- gen_post(fixed_parms=p0)   # fixed parameters can still be overriden by supplying them explicitly

print(lpost(p0))

ctrl <- list(fnscale=-1, maxit=2000)
opt_uncons <- optim(pstrt, lpost, control=ctrl)
puncons <- opt_uncons$par

pstrt2 <- puncons[-c(3)]
opt_nomob <- optim(pstrt2, lpost, control=ctrl)
pnomob <- opt_nomob$par
pnomob <- c(pnomob[1:2], zeta=0, pnomob[3:8])

## Looks like there's some local optimum stuff happening here.  Try another 
## optimization with zeta in effect
pstrt3 <- pnomob
opt_uncons2 <- optim(pstrt3, lpost, control=ctrl)
puncons2 <- opt_uncons2$par

## Results of the optimization calcs:
# puncons <-
# c(eta = -1.20312381497662, xi = 0.0669171609667788, zeta = 0.50413715543652, 
# D0 = 4.38140679701526, A0 = 0.711012863788961, b = 16.9753332901456, 
# I0 = 18.1657210691149, Ts = 4.18287942808935, mask_effect = -0.239518854477301
# )
# pnomob <-
# c(eta = -0.906737350533472, xi = 0.0392680308561262, zeta = 0, 
# D0 = 2.72167163219254, A0 = 1.09467383518768, b = 47.9638204551691, 
# I0 = 47.0287864607171, Ts = 3.18676012059561, mask_effect = -0.105449081847722
# )
# puncons2 <-
#   c(eta = -0.900792316608222, xi = 0.0340347134431018, zeta = -0.347512337099921, 
#     D0 = 2.61399709450331, A0 = 1.22579806178686, b = 46.4842688663428, 
#     I0 = 63.3532328749135, Ts = 2.85601231280208, mask_effect = -0.0975428701526071
#   )

pmat <- rbind(puncons, pnomob)
scen_names <- c('Mobility adjust', 'No mobility adjust')

proj <- plt_projections(pmat, scen_names, p0, usedate=TRUE, what='PopInfection') + scale_color_brewer(type='qual')
print(proj + ylab('Infected Population'))

modobs <- plt_modobs(pmat, scen_names, viscounties) + ggplot2::scale_color_brewer(type='qual')
print(modobs)

