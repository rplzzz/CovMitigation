library(metrosamp)
library(CovMitigation)
library(doParallel)


#pstrt <- c(T0 = 3.1, D0 = 4.9, A0 = 3.5, day_zero = 37.1, b = 3.1, I0=2.7)
pstrt <- c(T0 = 7.5, D0 = 4.5, A0 = 5.0, day_zero = 3.9, b = 8.7, I0=24.8)
scl <- c(0.1, 0.1, 0.1, 1.0, 1.0)

lpost <- gen_post()

run_mcmc <- function(tid, nsamp, outfilename, restartfile, nproc=8)
{
  if(is.na(nproc) || nproc==1) {
    foreach::registerDoSEQ()
  } else {
    registerDoParallel(cores=nproc)
  }
  
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
