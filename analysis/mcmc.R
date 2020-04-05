library(CovMitigation)
library(metrosamp)
library(doParallel)

registerDoParallel(8)

# Ts is not identifiable from this data, so fix it at its default value.
p0 <- c(T0=7, D0=4, A0=3, day_zero=60, b=20)
scl0 <- c(0.1, 0.1, 0.1, 1, 1)

lpost <- gen_post()

print(lpost(p0))

# Optimize with D0 and A0 fixed at their default values
pstrt <- c(T0=7, day_zero=50, b=5)
opt <- optim(pstrt, lpost, control=list(fnscale=-1))


set.seed(867-5309)


## Based on results of opt, above (caching these values to prevent having to 
## rerun it, since it takes about 45 minutes)
p1 <- c(T0=2.938418, D0=4, A0=3, day_zero=49.015829, b=10.304016)

ms1 <- metrosamp(lpost, p1, 100, 1, scl0)
ms2 <- metrosamp(lpost, ms1, 100,1, scl0/5)
ms3 <- metrosamp(lpost, ms2, 100, 1, scl0/2)
ms4 <- metrosamp(lpost, ms3, 100, 1, scl0/2)

ms5 <- metrosamp(lpost, ms4, 1000, 1, scl0/2)
ms6 <- metrosamp(lpost, ms5, 1000, 1, scl0)

