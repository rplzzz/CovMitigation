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
pstrt <- c(D0=7.1, A0=5, I0=10, Ts=4,  day_zero=10, b=10)
opt_fix_all_t0 <- optim(pstrt, lpost, control=list(fnscale=-1, maxit=1000))

pstrt_fix_gen_t0 <- c(T0_hi=7.5, opt_fix_all_t0$par)
opt_fix_gen_t0 <- optim(pstrt_fix_gen_t0, lpost, control=list(fnscale=-1, maxit=1000))

opt_uncons <- optim(c(T0=7.5, opt_fix_gen_t0$par), lpost, control=list(fnscale=-1, maxit=1000))

## Run an optimization with D0 fixed at something longer.
fplongd <- p0
fplongd['D0'] <- 10
pstrtlongd <- opt_uncons$par[-which(names(opt_uncons$par) == 'D0')]
lpost_longd <- gen_post(fixed_parms = fplongd)
opt_longd <- optim(pstrtlongd, lpost_longd, control=list(fnscale=-1, maxit=1000))

## run an optimization with b fixed at 100
fphib <- p0
fphib['b'] <- 50
pstrthib <- opt_uncons$par[-which(names(opt_uncons$par) == 'b')]
lpost_hib <- gen_post(fixed_parms = fphib)
opt_hib <- optim(pstrthib, lpost_hib, control=list(fnscale=-1, maxit=1000))


## Based on results of opt, above (caching these values to prevent having to 
## rerun it, since it takes kind of a long time.)
popt_t0_all <- c(T0_hi=7.5, T0=7.5, D0=1.8, A0=2.7, I0=19.0, Ts=1.9, day_zero=9.6, b=11.1)
popt_t0_gen <- c(T0_hi=6.3, T0=7.5, D0=1.9, A0=3.0, I0=4.2, Ts=1.9, day_zero=0.003, b=10.0)
popt_uncons <- c(T0_hi=6.3, T0=7.5, D0=1.9, A0=3.1, I0=4.3, Ts=1.9, day_zero=0.004, b=9.9)
popt_longd <- c(T0_hi=6.7, T0=8.4, D0=10, A0=3.6, I0=1.8, Ts=4.4, day_zero=0.002, b=12.9)
popt_hib <- c(T0_hi=7.0, T0=8.9, D0=2.4, A0=3.3, I0=1.8, Ts=2.2, day_zero=6.2, b=50)

pmat <- rbind(popt_uncons, popt_t0_all, popt_t0_gen, popt_longd, popt_hib)
scen_names <- 
  c( 'Unconstrained', 'Constrained td=7.5', 
   'Constrained td=7.5 (ex. NOVA)',
   'Constrained D0=10', 'Constrained b=50')
proj <- plt_projections(pmat, scen_names, usedate=TRUE)
print(proj + scale_color_brewer(type='qual'))

pltcounties <- c('FairfaxCounty', 'PrinceWilliamCounty', 'AlbemarleCounty', 'Charlottesvillecity')
modobs <- plt_modobs(pmat, scen_names, pltcounties)
print(modobs + ggplot2::facet_wrap(~locality, scales='free_y') + scale_color_brewer(type='qual'))

scl0 <- c(0.1, 0.1, 0.1, 1, 1, 1, 0.1)

plot_traces <- function(ms) {
  pltdata <- as.data.frame(ms$samples)
  pltdata$iter <- seq(1,nrow(pltdata))
  pltdata <- pivot_longer(pltdata, -iter, names_to = 'parm')
  ggplot(pltdata, aes(x=iter, y=value)) + geom_line() + facet_wrap(~parm, scales='free')
}

set.seed(867-5309)
ms1 <- metrosamp(lpost, popt_uncons, 100, 1, scl0)
ms2 <- metrosamp(lpost, ms1, 100,1, scl0/5)
ms3 <- metrosamp(lpost, ms2, 100, 1, scl0/5)
ms4 <- metrosamp(lpost, ms3, 100, 1, scl0/2)
ms5 <- metrosamp(lpost, ms4, 100, 1, scl0/2)
ms6 <- metrosamp(lpost, ms5, 100, 1, scl0/5)

ms7 <- metrosamp(lpost, ms6, 1000, 1, scl0/2)

