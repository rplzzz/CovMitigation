library(metrosamp)
library(CovMitigation)
library(doParallel)

registerDoParallel(40)

pstrt <- c(T0 = 2.938418, D0 = 4, A0 = 3, day_zero = 49.015829, b = 10.304016)

scl <- c(0.05, 0.05, 0.05, 0.50, 0.50)

lpost <- gen_post()

run_mcmc <- function(tid, nsamp, outfilename, restartfile)
{
  set.seed(867-5309 + tid)
  if(is.null(restartfile)) {
    mcs <- metrosamp(lpost, pstrt, nsamp, 1, scl)
  }
  else {
    warmup <- readRDS(restartfile)
    mcs <- metrosamp(lpost, warmup, nsamp,1)
  }
  saveRDS(mcs, outfilename, compress='xz')
}
