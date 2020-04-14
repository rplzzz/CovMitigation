library(metrosamp)
library(CovMitigation)
library(doParallel)


pstrt <- c(T0 = 5.0, D0 = 4.9, A0 = 3.5,  I0=3.3, Ts=3.6, day_zero = 6.2, b = 1.3)
scl <- c(0.05, 0.05, 0.05, 0.5, 0.5, 0.5, 0.05)

scl <- structure(c(0.025, -0.0073, -0.00033, 0.012, 0.0016, -0.31, 0.0079, 
-0.0073, 0.065, -0.00059, -0.003, 0.018, 0.031, -0.018, -0.00033, 
-0.00059, 0.028, -0.00065, -0.00012, -0.02, -0.00082, 0.012, 
-0.003, -0.00065, 0.18, 0.0018, 0.21, 0.0031, 0.0016, 0.018, 
-0.00012, 0.0018, 0.02, -0.055, -0.0046, -0.31, 0.031, -0.02, 
0.21, -0.055, 5.1, -0.068, 0.0079, -0.018, -0.00082, 0.0031, 
-0.0046, -0.068, 0.008), .Dim = c(7L, 7L), .Dimnames = list(c("T0", 
"D0", "A0", "I0", "Ts", "day_zero", "b"), c("T0", "D0", "A0", 
"I0", "Ts", "day_zero", "b")))

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
  
  if(!is.null(restartfile) && outfilename == restartfile) {
    file.rename(restartfile, namebackup(restartfile))
  }
  saveRDS(mcs, outfilename, compress='xz')
}
