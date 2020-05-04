library(metrosamp)
library(CovMitigation)
library(doParallel)

scl_warmup <-
  structure(c(1.441, 1.167, 10.03, -1.925, -0.8166, -0.4529, -0.2993, 
              5.888, 4.792, 1.167, 0.9665, 8.042, -2.045, -0.6536, -0.3639, 
              -0.2359, 4.591, 3.948, 10.03, 8.042, 100, -5.292, -4.696, -2.725, 
              -1.412, 60.55, 5.403, -1.925, -2.045, -5.292, 100, 4.204, 1.135, 
              1.263, 8.568, 37.21, -0.8166, -0.6536, -4.696, 4.204, 1.227, 
              0.3574, 0.5296, -2.47, -1.95, -0.4529, -0.3639, -2.725, 1.135, 
              0.3574, 0.1927, 0.1486, -1.449, -1.473, -0.2993, -0.2359, -1.412, 
              1.263, 0.5296, 0.1486, 0.2814, -0.641, -1.277, 5.888, 4.591, 
              60.55, 8.568, -2.47, -1.449, -0.641, 46.14, -4.437, 4.792, 3.948, 
              5.403, 37.21, -1.95, -1.473, -1.277, -4.437, 219.3), 
            .Dim = c(9L, 9L), 
            .Dimnames = list(c("T0_uhi", "T0_hi", "T0_lo", "T0_ulo", 
                               "D0", "A0", "Ts", "I0", "b"), 
                             c("T0_uhi", "T0_hi", "T0_lo", "T0_ulo", 
                               "D0", "A0", "Ts", "I0", "b")))

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
lpost <- gen_post(hparms=hyper_parms, fixed_parms = c(day_zero=30))

## We generally run 16 independent chains.  The only time this exact number matters
## is when searching for an initial guess during the warmup runs.
nchain <- 16

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
    
    ## Search for a starting parameter set using quasi-random draws from the priors.
    ## Occasionally we will draw a set that doesn't produce a finite log-posterior, 
    ## so keep on searching until we find one that is valid.
    irow <- 1 + tid    
    repeat {
      qrvals <- randtoolbox::sobol(irow, dim=ndim)
      pstrt <- qprior(qrvals[irow,], hyper_parms)
      if(pstrt['day_zero'] <= 0) {
        ## Illegal parameter value -- shift to > 0
        pstrt['day_zero'] <- 0.01
      }
      lpval <- lpost(pstrt)
      sstrt <- paste(names(pstrt), " :\t", signif(pstrt,3), collapse='\n')
      message('Proposed start:\n', sstrt, '\nlpost = ', lpval, '\n')
      if(is.finite(lpval)) {
        break
      }
      irow <- irow + nchain
    }
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
