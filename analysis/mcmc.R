library(CovMitigation)
library(tidyr)
library(dplyr)
library(ggplot2)
library(metrosamp)
library(doParallel)

registerDoParallel(8)

viscounties <- c('Harrisonburgcity', 'AugustaCounty',   ## ultra-high
                 'FairfaxCounty', 'HenricoCounty',      ## high
                 'AlbemarleCounty', 'Charlottesvillecity', ## low
                 'JamesCityCounty', 'AmherstCounty'
)


p0 <- c(T0_uhi = 5, T0_hi = 7.5, T0_lo=14, T0_ulo=25, D0=7, A0=5, day_zero=60, b=20, I0=10, Ts=4)

lpost <- gen_post(fixed_parms=p0)   # fixed parameters can still be overriden by supplying them explicitly

print(lpost(p0))

# Optimize with T0 fixed at its default value
ctrl <- list(fnscale=-1, maxit=2000)
pstrt <- c(D0=7.1, A0=5, I0=10, Ts=4,  day_zero=10, b=10)
opt_fix_all_t0 <- optim(pstrt, lpost, control=ctrl)
pfixallt0 <- c(p0[1:4], opt_fix_all_t0$par)

opt_uncons <- optim(pfixallt0, lpost, control=ctrl)
puncons <- opt_uncons$par

## Run an optimization with D0 fixed at something longer.
fplongd <- puncons
fplongd['D0'] <- 10
pstrtlongd <- opt_uncons$par[-which(names(opt_uncons$par) == 'D0')]
lpost_longd <- gen_post(fixed_parms = fplongd)
opt_longd <- optim(pstrtlongd, lpost_longd, control=ctrl)

## run an optimization with b fixed at 100
fphib <- puncons
fphib['b'] <- 50
pstrthib <- opt_uncons$par[-which(names(opt_uncons$par) == 'b')]
lpost_hib <- gen_post(fixed_parms = fphib)
opt_hib <- optim(pstrthib, lpost_hib, control=ctrl)

## run an optimized "pause" scenario, with tau fixed at 150 days
fppause <- puncons
fppause['T0_uhi'] <- 150
fppause['T0_hi'] <- 150
fppause['T0_lo'] <- 150
fppause['T0_ulo'] <- 150
pstrtpause <- opt_uncons$par[-which(grepl('T0', names(opt_uncons$par)))]
pstrtpause['day_zero'] <- 60
pstrtpause['b'] <- 100
lpost_pause <- gen_post(fixed_parms = fppause)
opt_pause <- optim(pstrtpause, lpost_pause, control=ctrl)

## Another pause scenario, but in this one only the low and ultra-low rates
## are fixed
fppause2 <- puncons
fppause2['T0_lo'] <- 150
fppause2['T0_ulo'] <- 150
pstrtpause2 <- opt_uncons$par[-which(grepl('lo', names(opt_uncons$par)))]
pstrtpause2['day_zero'] <- 60
pstrtpause2['b'] <- 100
lpost_pause2 <- gen_post(fixed_parms = fppause2)
opt_pause2 <- optim(pstrtpause2, lpost_pause2, control=ctrl)


## Based on results of opt, above (caching these values to prevent having to 
## rerun it, since it takes kind of a long time.)
popt_t0_all <- c(T0_uhi=5, T0_hi=7.5, T0_lo=15, T0_ulo=25, D0=2.0, A0=0.68, I0=4.26, Ts=3.7, day_zero=15.2, b=24.0)
popt_uncons <- c(T0_uhi=11.5, T0_hi=9.1, T0_lo=15.4, T0_ulo=28.3, D0=1.9, A0=1.4, I0=9.7, Ts=2.5, day_zero=0.002, b=27.2)
popt_longd <- c(T0_uhi=10.9, T0_hi=9.6, T0_lo=14.1, T0_ulo=23.2, D0=10, A0=2.8, I0=3.0, Ts=5.1, day_zero=9e-5, b=27.2)
popt_hib <- c(T0_uhi=11.8, T0_hi=9.7, T0_lo=16.4, T0_ulo=30.1, D0=1.3, A0=1.4, I0=10.8, Ts=1.6, day_zero=0.003, b=50)
popt_pause <- c(T0_uhi=150, T0_hi=150, T0_lo=150, T0_ulo=150, D0=2.0, A0=1.0, I0=214, Ts=2.4, day_zero=2.3, b=354)
popt_pause2 <- c(T0_uhi=10.9, T0_hi=8.5, T0_lo=150, T0_ulo=150, D0=7.2, A0=2.9, I0=25.5, Ts=4.6, day_zero=53.0, b=96.0)

#popt_uncons2 <- c(T0_hi=6.2, T0=9.5, D0=6.7, A0=3.8, I0=11.4, Ts=2.6, day_zero=43.1, b=91.1)

pmat <- rbind(popt_uncons, popt_t0_all, popt_longd, popt_hib, popt_pause,
              popt_pause2)
#, popt_uncons2)
scen_names <- 
  c( 'Unconstrained', 'Constrained growth rates', 
     'Constrained D0=10', 'Constrained b=50',
     'Pause scenario', 'Pause scenario (ex. hi and uhi)')
#,'Alternate Unconstrained')

lfcmp <- sapply(seq(1,nrow(pmat)), 
                function(i) {
                  lpost(pmat[i,])
                })
ii <- order(lfcmp)
for(i in ii) {
  cat(paste(scen_names[i], ':\t', lfcmp[i], '\n'))
}

proj <- plt_projections(pmat, scen_names, usedate=TRUE)
print(proj + ggplot2::scale_color_brewer(type='qual') + 
        ggplot2::xlim(as.Date(c('2020-03-20', '2020-06-01'))))

modobs <- plt_modobs(pmat, scen_names, viscounties)
print(modobs + ggplot2::scale_color_brewer(type='qual'))

scl0 <- c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 1, 1, 1, 0.1)

plot_traces <- function(ms) {
  pltdata <- as.data.frame(ms$samples)
  pltdata$iter <- seq(1,nrow(pltdata))
  pltdata <- pivot_longer(pltdata, -iter, names_to = 'parm')
  ggplot(pltdata, aes(x=iter, y=value)) + geom_line() + facet_wrap(~parm, scales='free')
}

set.seed(867-5309)
ms1 <- metrosamp(lpost, popt_hib, 100, 1, scl0)
ms2 <- metrosamp(lpost, ms1, 100,1, scl0/2)
ms3 <- metrosamp(lpost, ms2, 100, 1, scl0/5)
ms4 <- metrosamp(lpost, ms3, 100, 1, scl0/5)
ms5 <- metrosamp(lpost, ms4, 100, 1, scl0/5)

ms6 <- metrosamp(lpost, ms5, 1000, 1, scl0/5)
ms7 <- metrosamp(lpost, ms6, 1000, 1, scl0/10)


