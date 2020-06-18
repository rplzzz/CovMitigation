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

day0fix <- 30
p0 <- c(T0_uhi = 5, T0_hi = 7.5, T0_lo=14, T0_ulo=25, D0=7, A0=5, day_zero=day0fix, b=20, I0=10, Ts=4)

lpost <- gen_post(fixed_parms=p0)   # fixed parameters can still be overriden by supplying them explicitly

print(lpost(p0))

# Optimize with T0 fixed at its default value
ctrl <- list(fnscale=-1, maxit=2000)
pstrt <- c(D0=7.1, A0=5, I0=10, Ts=4, b=10)
opt_fix_all_t0 <- optim(pstrt, lpost, control=ctrl)
pfixallt0 <- c(p0[1:4], opt_fix_all_t0$par)

opt_uncons <- optim(pfixallt0, lpost, control=ctrl)
puncons <- opt_uncons$par

## Run an optimization with D0 fixed at something longer.
fplongd <- puncons
fplongd['D0'] <- 7
pstrtlongd <- opt_uncons$par[-which(names(opt_uncons$par) == 'D0')]
lpost_longd <- gen_post(fixed_parms = c(fplongd, day_zero=day0fix))
opt_longd <- optim(pstrtlongd, lpost_longd, control=ctrl)

## run an optimization with b fixed at 100
fphib <- puncons
fphib['b'] <- 50
pstrthib <- opt_uncons$par[-which(names(opt_uncons$par) == 'b')]
lpost_hib <- gen_post(fixed_parms = c(fphib, day_zero=day0fix))
opt_hib <- optim(pstrthib, lpost_hib, control=ctrl)

## run an optimized "pause" scenario, with tau fixed at 150 days
fppause <- puncons
fppause['T0_uhi'] <- 150
fppause['T0_hi'] <- 150
fppause['T0_lo'] <- 150
fppause['T0_ulo'] <- 150
pstrtpause <- opt_uncons$par[-which(grepl('T0', names(opt_uncons$par)))]
pstrtpause['b'] <- 100
lpost_pause <- gen_post(fixed_parms = c(fppause, day_zero=day0fix))
opt_pause <- optim(pstrtpause, lpost_pause, control=ctrl)

## Another pause scenario, but in this one only the low and ultra-low rates
## are fixed
fppause2 <- puncons
fppause2['T0_lo'] <- 150
fppause2['T0_ulo'] <- 150
pstrtpause2 <- opt_uncons$par[-which(grepl('lo', names(opt_uncons$par)))]
pstrtpause2['b'] <- 100
lpost_pause2 <- gen_post(fixed_parms = c(fppause2, day_zero=day0fix))
opt_pause2 <- optim(pstrtpause2, lpost_pause2, control=ctrl)


## Based on results of opt, above (caching these values to prevent having to 
## rerun it, since it takes kind of a long time.)
popt_t0_all <- c(T0_uhi=5, T0_hi=7.5, T0_lo=15, T0_ulo=25, D0=2.9, A0=0.78, I0=7.1, Ts=4.9, b=11.6)
popt_uncons <- c(T0_uhi=10.5, T0_hi=8.0, T0_lo=16.1, T0_ulo=4.1, D0=2.4, A0=1.1, I0=13.2, Ts=2.7, b=16.2)
popt_longd <- c(T0_uhi=13.6, T0_hi=10.2, T0_lo=18.3, T0_ulo=3.9, D0=10, A0=1.9, I0=35.7, Ts=5.5, b=8.0)
popt_hib <- c(T0_uhi=10.4, T0_hi=9.1, T0_lo=17.8, T0_ulo=5.0, D0=1.4, A0=0.96, I0=13.8, Ts=1.6, b=50)
popt_pause <- c(T0_uhi=150, T0_hi=150, T0_lo=150, T0_ulo=150, D0=5.7, A0=2.7, I0=66.1, Ts=5.4, b=176)
popt_pause2 <- c(T0_uhi=33.3, T0_hi=19.6, T0_lo=150, T0_ulo=150, D0=3.9, A0=2.2, I0=38.3, Ts=3.8, b=76)

pmat <- rbind(popt_uncons, popt_t0_all, popt_longd, popt_hib, popt_pause,
              popt_pause2)

scen_names <- 
  c( 'Unconstrained', 'Constrained growth rates', 
     'Constrained D0=10', 'Constrained b=50',
     'Pause scenario', 'Pause scenario (ex. hi and uhi)')

lfcmp <- sapply(seq(1,nrow(pmat)), 
                function(i) {
                  lpost(pmat[i,])
                })
ii <- order(lfcmp)
for(i in ii) {
  cat(paste(scen_names[i], ':\t', lfcmp[i], '\n'))
}

proj <- plt_projections(pmat, scen_names, p0, usedate=TRUE)
print(proj + ggplot2::scale_color_brewer(type='qual') + 
        ggplot2::xlim(as.Date(c('2020-04-15', '2020-07-31'))))

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


