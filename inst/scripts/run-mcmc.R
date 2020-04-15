library(metrosamp)
library(CovMitigation)
library(doParallel)


pstrt <- c(T0 = 5.0, D0 = 4.9, A0 = 3.5,  I0=3.3, Ts=3.6, day_zero = 6.2, b = 1.3)

pmat <- 
  structure(c(6.3, 7.5, 6.3, 6.7, 7, 7.5, 7.5, 7.5, 8.4, 8.9, 1.9, 
              1.8, 1.9, 10, 2.4, 3.1, 2.7, 3, 3.6, 3.3, 4.3, 19, 4.2, 1.8, 
              1.8, 1.9, 1.9, 1.9, 4.4, 2.2, 0.004, 9.6, 0.003, 0.002, 6.2, 
              9.9, 11.1, 10, 12.9, 50), .Dim = c(5L, 8L), .Dimnames = list(
                c("popt_uncons", "popt_t0_all", "popt_t0_gen", "popt_longd", 
                  "popt_hib"), c("T0_hi", "T0", "D0", "A0", "I0", "Ts", "day_zero", 
                                 "b")))

scl <-
  structure(c(0.025, 0, 0, 0, 0, 0, 0, 0, 0, 0.025, -0.0073, -0.00033, 
              0.012, 0.0016, -0.31, 0.0079, 0, -0.0073, 0.065, -0.00059, -0.003, 
              0.018, 0.031, -0.018, 0, -0.00033, -0.00059, 0.028, -0.00065, 
              -0.00012, -0.02, -0.00082, 0, 0.012, -0.003, -0.00065, 0.18, 
              0.0018, 0.21, 0.0031, 0, 0.0016, 0.018, -0.00012, 0.0018, 0.02, 
              -0.055, -0.0046, 0, -0.31, 0.031, -0.02, 0.21, -0.055, 5.1, -0.068, 
              0, 0.0079, -0.018, -0.00082, 0.0031, -0.0046, -0.068, 0.008), 
            .Dim = c(8L, 8L), 
            .Dimnames = list(c("T0_hi", "T0", "D0", "A0", "I0", "Ts", "day_zero", "b"), 
                             c("T0_hi", "T0", "D0", "A0", "I0", "Ts", "day_zero", "b")))


lpost <- gen_post()

run_mcmc <- function(tid, nsamp, outfilename, restartfile=NULL, 
                     nproc=8, nwarmup = 1000, usescl=TRUE)
{
  if(is.na(nproc) || nproc==1) {
    foreach::registerDoSEQ()
  } else {
    registerDoParallel(cores=nproc)
  }
  
  set.seed(867-5309 + tid)
  if(is.null(restartfile)) {
    ## Distribute the samplers amongst the optimization results.
    irow <- 1 + tid%%nrow(pmat)    
    pstrt <- pmat[irow,]
    
    mcwarmup <- metrosamp(lpost, pstrt, nwarmup, 1, scl)
    mcs <- metrosamp(lpost, mcwarmup, nsamp, 1)
  }
  else {
    mcwarmup <- readRDS(restartfile)
    if(usescl) {
      mcs <- metrosamp(lpost, mcwarmup, nsamp, 1, scl)
    } else {
      mcs <- metrosamp(lpost, mcwarmup, nsamp,1)
    }
  }
  
  if(!is.null(restartfile) && outfilename == restartfile) {
    file.rename(restartfile, namebackup(restartfile))
  }
  saveRDS(mcs, outfilename, compress='xz')
  
  cat('\nFIN.\n')
}
