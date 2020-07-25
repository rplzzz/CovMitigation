library(CovMitigation)
library(metrosamp)
library(doParallel)

lpost <- gen_post()
nchain <- 16
p0 <- gen_parm_vec()                  # Parameter vector for all counties
ndim <- length(p0)

mcmc_init <- function(tid, outfilename, nproc=8)
{
  if(is.na(nproc) || nproc==1) {
    foreach::registerDoSEQ()
  } else {
    registerDoParallel(cores=nproc)
  }
  
  set.seed(867-5309 + tid)
  
  ## Search for a starting parameter set using quasi-random draws from the priors.
  ## Occasionally we will draw a set that doesn't produce a finite log-posterior, 
  ## so keep on searching until we find one that is valid.
  irow <- 1 + tid    
  repeat {
    qrvals <- randtoolbox::sobol(irow, dim=ndim)
    pv <- qrvals[irow,]
    names(pv) <- names(p0)
    pstrt <- qprior(pv)
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
  opt <- optim(pstrt, lpost, control=list(fnscale=-1, maxit=5000))

  ## Convert these parameter values to a metrosamp structure by doing a single step
  scl <- rep(0.1, length(opt$par))
  mcs <- metrosamp(lpost, opt$par, 1, 1, scl, debug=TRUE)

  plast <- mcs$plast
  slast <- paste(names(plast), ' :\t', signif(plast,3), collapse='\n')
  lplast <- lpost(plast)
  message('Final optimized value:\n', slast, '\nlpost = ', lplast, '\n')
  
  saveRDS(mcs, outfilename, compress='xz')
  
  cat('\nFIN.\n')
}
