library(metrosamp)
library(CovMitigation)
library(doParallel)

registerDoParallel(20)

pstrt <-
c(T0 = 3.10253101176298, D0 = 8.09589183577002, A0 = 4.39695340708852, 
day_zero = 53.5297768069546, b = 75.1143457008878)

scl <- structure(c(0.25, 0, 0, -0.75, 0.75, 0, 0.25, 0, 0, 0, 0, 0, 
0.25, 0, 0, -0.75, 0, 0, 6.25, 0, 0.75, 0, 0, 0, 6.25), .Dim = c(5L, 
5L))

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
