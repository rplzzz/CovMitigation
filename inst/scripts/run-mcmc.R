library(metrosamp)
library(CovMitigation)
library(doParallel)

scl_warmup <- c(0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.1, 0.1, 0.1, 0.01)

scl_main <-
  structure(c(0.0196274353762327, 0.0296350461689918, 0, 0, -0.0101040739504531, 
              0.00433204814709052, -0.228175572395601, 0.0920387639847004, 
              0.0296350461689918, 0.057437930699364, -0.0148183235736357, 0, 
              0.0132642755736333, 0, -0.265366089216169, 0.183813393331308, 
              0, -0.0148183235736357, 0.29739627160971, 0, 0, 0.0351280178739933, 
              0, -1.81716470695934, 0, 0, 0, 0.0624104813033857, 0, 0, 0, 0, 
              -0.0101040739504531, 0.0132642755736333, 0, 0, 0.152675234085546, 
              0, 0.631655044192764, 0, 0.00433204814709052, 0, 0.0351280178739934, 
              0, 0, 0.0333033207758625, -0.068196234565824, -0.391281938547136, 
              -0.228175572395601, -0.265366089216169, 0, 0, 0.631655044192764, 
              -0.068196234565824, 4.78644609255565, 0, 0.0920387639847004, 
              0.183813393331308, -1.81716470695934, 0, 0, -0.391281938547136, 
              0, 22.2799178084844), 
            .Dim = c(8L, 8L), 
            .Dimnames = list(c("T0_hi","T0", "D0", "A0", "I0", "Ts", "day_zero", "b"), 
                             c("T0_hi", "T0", "D0", "A0", "I0", "Ts", "day_zero", "b")))

hyper_parms <- list(nhosp_weight=50)
lpost <- gen_post(hparms=hyper_parms)

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
    ## This calculation is a warmup.  Start with quasi-random samples from the
    ## prior distribution.
    if(is.matrix(scl_warmup)) {
      ndim <- ncol(scl_warmup)
    }
    else {
      ndim <- length(scl_warmup)
    }
    qrvals <- randtoolbox::sobol(tid+1, dim=ndim)
    irow <- 1 + tid    
    pstrt <- qprior(qrvals[irow,], hyper_parms)
    
    ## Now do a local optimization to find a "pretty good" start value.  It's 
    ## ok if this doesn't make it all the way to convergence.
    opt <- optim(pstrt, lpost, control=list(fnscale=-1))
    
    mcs <- metrosamp(lpost, opt$par, nsamp, 1, scl_warmup)
  }
  else {
    ## Production run, either starting from a warmup run or continuing a previous batch.
    mcwarmup <- readRDS(restartfile)
    if(usescl) {
      mcs <- metrosamp(lpost, mcwarmup, nsamp, 1, scl_main)
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
