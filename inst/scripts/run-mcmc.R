library(metrosamp)
library(CovMitigation)
library(doParallel)


pstrt <- c(T0 = 2.7, D0 = 3.6, A0 = 3.6, day_zero = 39.2, b = 1.1, I0=4.1, Ts=3.6)
scl <- c(0.05, 0.05, 0.05, 0.5, 0.5, 0.5, 0.05)

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
  
  if(outfilename == restartfile) {
    file.rename(restartfile, namebackup(restartfile))
  }
  saveRDS(mcs, outfilename, compress='xz')
}
