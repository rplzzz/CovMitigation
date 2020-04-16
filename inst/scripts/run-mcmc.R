library(metrosamp)
library(CovMitigation)
library(doParallel)


pmat <-
  structure(c(5.9, 7.5, 6.2, 6.8, 7.1, 150, 7.4, 6.2, 7.1, 7.5, 
              7.5, 8.8, 9.4, 150, 186, 9.5, 2.1, 1.9, 2.4, 10, 2.3, 9.8, 9.3, 
              6.7, 3.2, 2.8, 3.1, 3.5, 2.8, 2.6, 3.1, 3.8, 4.5, 19.6, 7.2, 
              3, 2.3, 149.7, 64.4, 11.4, 1.9, 2, 2.1, 4.2, 2.2, 5.3, 2.8, 2.6, 
              6.1, 11, 7, 4.1, 6.3, 43, 56.9, 43.1, 9.4, 11.2, 8.7, 13.8, 50, 
              96.8, 153.8, 91.1), .Dim = c(8L, 8L), 
            .Dimnames = list(c("popt_uncons", 
                               "popt_t0_all", "popt_t0_gen", "popt_longd", "popt_hib", "popt_pause", 
                               "popt_pause_exnova", "popt_uncons2"), 
                             c("T0_hi", "T0", "D0", 
                               "A0", "I0", "Ts", "day_zero", "b")))

scl <-c(0.05, 0.05, 0.05, 0.05, 0.5, 0.5, 0.5, 0.05)

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
