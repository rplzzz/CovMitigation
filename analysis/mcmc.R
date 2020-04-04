library(CovMitigation)
library(metrosamp)

# Ts is not identifiable from this data, so fix it at its default value.
p0 <- c(T0=7, D0=4, A0=3, day_zero=60, b=20)
scl0 <- c(0.1, 0.1, 0.1, 1, 1)

lpost <- gen_post()

print(lpost(p0))

set.seed(867-5309)

#opt <- optim(p0, lpost, control=list(fnscale=-1))

## Based on results of opt, above (caching these values to prevent having to 
## rerun it, since it takes about 45 minutes)
p1 <- c(T0=5.48, D0=0.08, A0=0.3, day_zero=25.2, b=46.3)

ms1 <- metrosamp(lpost, p1, 100, 1, scl0)
ms2 <- metrosamp(lpost, ms1, 100,1)
ms3 <- metrosamp(lpost, ms2, 100, 1)
ms4 <- metrosamp(lpost, ms3, 100, 1, scl0/5)

scl_mat0 <- cor2cov(cor(ms4$samples), scl0)
ms5 <- metrosamp(lpost, ms4, 100, 1, scl_mat0)
ms6 <- metrosamp(lpost, ms5, 100, 1)
