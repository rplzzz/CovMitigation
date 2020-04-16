library(CovMitigation)
library(tidyr)
library(dplyr)
library(ggplot2)
library(metrosamp)
library(doParallel)

registerDoParallel(8)

# Ts is not identifiable from this data, so fix it at its default value.
p0 <- c(T0_hi = 7.5, T0=7.5, D0=7, A0=5, day_zero=60, b=20, I0=10, Ts=4)
modp0 <- c(T0=7.5, D0=7, A0=5, I0=10, Ts=4)

lpost <- gen_post(fixed_parms=p0)   # fixed parameters can still be overriden by supplying them explicitly

print(lpost(p0))

# Optimize with T0 fixed at its default value
ctrl <- list(fnscale=-1, maxit=2000)
pstrt <- c(D0=7.1, A0=5, I0=10, Ts=4,  day_zero=10, b=10)
opt_fix_all_t0 <- optim(pstrt, lpost, control=ctrl)

pstrt_fix_gen_t0 <- c(T0_hi=7.5, opt_fix_all_t0$par)
opt_fix_gen_t0 <- optim(pstrt_fix_gen_t0, lpost, control=ctrl)

opt_uncons <- optim(c(T0=7.5, opt_fix_gen_t0$par), lpost, control=ctrl)

## Run an optimization with D0 fixed at something longer.
fplongd <- p0
fplongd['D0'] <- 10
pstrtlongd <- opt_uncons$par[-which(names(opt_uncons$par) == 'D0')]
lpost_longd <- gen_post(fixed_parms = fplongd)
opt_longd <- optim(pstrtlongd, lpost_longd, control=ctrl)

## run an optimization with b fixed at 100
fphib <- p0
fphib['b'] <- 50
pstrthib <- opt_uncons$par[-which(names(opt_uncons$par) == 'b')]
lpost_hib <- gen_post(fixed_parms = fphib)
opt_hib <- optim(pstrthib, lpost_hib, control=ctrl)

## run an optimized "pause" scenario, with tau fixed at 150 days
fppause <- p0
fppause['T0_hi'] <- 150
fppause['T0'] <- 150
pstrtpause <- opt_uncons$par[-which(names(opt_uncons$par) %in% c('T0_hi','T0'))]
pstrtpause['day_zero'] <- 60
pstrtpause['b'] <- 100
lpost_pause <- gen_post(fixed_parms = fppause)
opt_pause <- optim(pstrtpause, lpost_pause, control=ctrl)

## Found these starting parameters through experimentation.
pstrtpause_exnova <- c(T0_hi = 5, T0 = 150, D0 = 13.9, A0 = 3, I0 = 93, Ts = 4.2, 
                       day_zero = 50.2, b = 158.9)
opt_pause_exnova <- optim(pstrtpause_exnova, lpost_pause, control=ctrl)

## Based on results of opt, above (caching these values to prevent having to 
## rerun it, since it takes kind of a long time.)
popt_t0_all <- c(T0_hi=7.5, T0=7.5, D0=1.9, A0=2.8, I0=19.6, Ts=2.0, day_zero=11.0, b=11.2)
popt_t0_gen <- c(T0_hi=6.2, T0=7.5, D0=2.4, A0=3.1, I0=7.2, Ts=2.1, day_zero=7.0, b=8.7)
popt_uncons <- c(T0_hi=5.9, T0=7.1, D0=2.1, A0=3.2, I0=4.5, Ts=1.9, day_zero=6.1, b=9.4)
popt_longd <- c(T0_hi=6.8, T0=8.8, D0=10, A0=3.5, I0=3.0, Ts=4.2, day_zero=4.1, b=13.8)
popt_hib <- c(T0_hi=7.1, T0=9.4, D0=2.3, A0=2.8, I0=2.3, Ts=2.2, day_zero=6.3, b=50)
popt_pause <- c(T0_hi=150, T0=150, D0=9.8, A0=2.6, I0=149.7, Ts=5.3, day_zero=43.0, b=96.8)
popt_pause_exnova <- c(T0_hi=7.4, T0=186, D0=9.3, A0=3.1, I0=64.4, Ts=2.8, day_zero=56.9, b=153.8)
popt_uncons2 <- c(T0_hi=6.2, T0=9.5, D0=6.7, A0=3.8, I0=11.4, Ts=2.6, day_zero=43.1, b=91.1)

pmat <- rbind(popt_uncons, popt_t0_all, popt_t0_gen, popt_longd, popt_hib, popt_pause,
              popt_pause_exnova, popt_uncons2)
scen_names <- 
  c( 'Unconstrained', 'Constrained td=7.5', 
   'Constrained td=7.5 (ex. NOVA)',
   'Constrained D0=10', 'Constrained b=50',
   'Pause scenario', 'Pause scenario (ex. NOVA)',
   'Alternate Unconstrained')

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

pltcounties <- c('FairfaxCounty', 'PrinceWilliamCounty', 'AlbemarleCounty', 'Charlottesvillecity')
modobs <- plt_modobs(pmat, scen_names, pltcounties)
print(modobs + ggplot2::facet_wrap(~locality, scales='free_y') + ggplot2::scale_color_brewer(type='qual'))

scl0 <- c(0.1, 0.1, 0.1, 0.1, 1, 1, 1, 0.1)

plot_traces <- function(ms) {
  pltdata <- as.data.frame(ms$samples)
  pltdata$iter <- seq(1,nrow(pltdata))
  pltdata <- pivot_longer(pltdata, -iter, names_to = 'parm')
  ggplot(pltdata, aes(x=iter, y=value)) + geom_line() + facet_wrap(~parm, scales='free')
}

set.seed(867-5309)
ms1 <- metrosamp(lpost, popt_uncons2, 100, 1, scl0)
ms2 <- metrosamp(lpost, ms1, 100,1, scl0*2)
ms3 <- metrosamp(lpost, ms2, 100, 1, scl0)
ms4 <- metrosamp(lpost, ms3, 100, 1, scl0/5)
ms5 <- metrosamp(lpost, ms4, 100, 1, scl0/2)

ms6 <- metrosamp(lpost, ms5, 1000, 1, scl0/2)

