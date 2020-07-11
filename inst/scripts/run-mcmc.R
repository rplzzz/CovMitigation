library(metrosamp)
library(CovMitigation)
library(doParallel)

scl_warmup <-
c(eta = 0.1, xi = 0.1, zeta = 0.1, D0 = 0.5, A0 = 0.5, b = 5, 
I0 = 1, Ts = 0.5, mask_effect = 0.1)

scl_main <- scl_warmup
  
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
      if('day_zero' %in% names(pstrt) && pstrt['day_zero'] <= 0) {
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
    
    mcs <- metrosamp(lpost, opt$par, nsamp, 1, scl_warmup, debug=TRUE)
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
