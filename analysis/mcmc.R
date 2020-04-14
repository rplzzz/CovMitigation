library(CovMitigation)
library(tidyr)
library(dplyr)
library(ggplot2)
library(metrosamp)
library(doParallel)

registerDoParallel(8)

# Ts is not identifiable from this data, so fix it at its default value.
p0 <- c(T0=7.5, D0=7, A0=5, day_zero=60, b=20, I0=10, Ts=4)
modp0 <- c(T0=7.5, D0=7, A0=5, I0=10, Ts=4)

lpost <- gen_post(fixed_parms=p0)   # fixed parameters can still be overriden by supplying them explicitly

print(lpost(p0))

# Optimize with T0 fixed at its default value
pstrt <- c(D0=7.1, A0=5, I0=10, Ts=4,  day_zero=10, b=10)
opt <- optim(pstrt, lpost, control=list(fnscale=-1, maxit=1000))

opt_uncons <- optim(c(T0=7.5, opt$par), lpost, control=list(fnscale=-1, maxit=1000))

## Run an optimization with D0 fixed at something longer.
fplongd <- p0
fplongd['D0'] <- 10
pstrtlongd <- opt_uncons$par[-which(names(opt_uncons$par) == 'D0')]
lpost_longd <- gen_post(fixed_parms = fplongd)
opt_longd <- optim(pstrtlongd, lpost_longd, control=list(fnscale=-1, maxit=1000))

## run an optimization with b fixed at 100
fphib <- p0
fphib['b'] <- 100
pstrthib <- opt_uncons$par[-which(names(opt_uncons$par) == 'b')]
lpost_hib <- gen_post(fixed_parms = fphib)
opt_hib <- optim(pstrthib, lpost_hib, control=list(fnscale=-1, maxit=1000))


## Based on results of opt, above (caching these values to prevent having to 
## rerun it, since it takes kind of a long time.)
popt <- c(T0=7.5, D0=4.0, A0=3.0, I0=21.1, Ts=3.5, day_zero=0.0006, b=3.0 )
popt_uncons <- c(T0=4.8, D0=4.8, A0=3.4, I0=3.3, Ts=3.4, day_zero=9.9, b=1.3)
popt_longd <- c(T0=4.5, D0=10, A0=3.5, I0=3.3, Ts=4.5, day_zero=10.9, b=0.5)
popt_hib <- c(T0=5.0, D0=6.3, A0=3.5, I0=2.1, Ts=3.0, day_zero=46, b=100)

pmat <- rbind(popt, popt_longd, popt_uncons, popt_hib)
plt_projections(pmat, c('Constrained td=7.5', 'Constrained D0=10', 'Unconstrained', 'Constrained b=100'), usedate=TRUE)

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

