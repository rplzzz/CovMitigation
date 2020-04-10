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
pstrt <- c(D0=7, day_zero=60, b=20, I0=10)
opt <- optim(pstrt, lpost, control=list(fnscale=-1))


set.seed(867-5309)


## Based on results of opt, above (caching these values to prevent having to 
## rerun it, since it takes about 45 minutes)
popt <- c(T0=7.5, D0=4.5, A0=5, day_zero=3.9, b=8.7, Ts=4, I0=24.8)

pmat <- matrix(c(p0,popt), nrow=2, byrow = TRUE)
plt_projections(pmat, c('Base Scenario', 'Opt Scenario'), usedate=TRUE)

p1 <- c(T0=7.5, D0=4.5, A0=5, day_zero=3.9, b=8.7, I0=24.8)
scl0 <- c(0.1, 0.1, 0.1, 1, 1, 1)

plot_traces <- function(ms) {
  pltdata <- as.data.frame(ms$samples)
  pltdata$iter <- seq(1,nrow(pltdata))
  pltdata <- pivot_longer(pltdata, -iter, names_to = 'parm')
  ggplot(pltdata, aes(x=iter, y=value)) + geom_line() + facet_wrap(~parm, scales='free')
}

ms1 <- metrosamp(lpost, p1, 100, 1, scl0)
ms2 <- metrosamp(lpost, ms1, 100,1, scl0*2)
ms3 <- metrosamp(lpost, ms2, 100, 1, scl0*5)
ms4 <- metrosamp(lpost, ms3, 100, 1, scl0*2)
ms5 <- metrosamp(lpost, ms4, 100, 1, scl0)

ms6 <- metrosamp(lpost, ms5, 1000, 1, scl0)
scl1 <- scl0
scl1[4] <- 5
scl1[5] <- 5
scl1[6] <- 2
ms7 <- metrosamp(lpost, ms6, 1000, 1, scl1)

