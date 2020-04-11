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

# Optimize with D0 and A0 fixed at their default values
pstrt <- c(D0=7, A0=5, day_zero=10, b=10, I0=10)
opt <- optim(pstrt, lpost, control=list(fnscale=-1, maxit=1000))

opt_uncons <- optim(c(T0=7.5, opt$par), lpost, control=list(fnscale=-1, maxit=1000))

## Based on results of opt, above (caching these values to prevent having to 
## rerun it, since it takes about 45 minutes)
popt <- c(T0=7.5, D0=3.4, A0=2.9, day_zero=5.9, b=4.13, I0=22.1, Ts=4)
popt_uncons <- c(T0=2.7, D0=3.6, A0=3.6, day_zero=39.2, b=1.1, I0=4.1, Ts=4)

pmat <- rbind(popt, popt_uncons)
plt_projections(pmat, c('Constrained td=7.5', 'Unconstrained'), usedate=TRUE)

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

ms6 <- metrosamp(lpost, ms5, 1000, 1)

