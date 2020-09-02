library(CovMitigation)
library(metrosamp)


ptrans <- function(m, icol, newcolname)
{
  xcdf <- ecdf(m[ ,icol])
  m[ ,icol] <- xcdf(m[ ,icol])
  colnames(m)[icol] <- newcolname

  m
}

read1ensemble <- function(filename)
{
  stopifnot(file.exists(filename))
  ms <- readRDS(filename)
  samps <- getsamples(ms)

  ## perform probability transform on betacol.  This allows us to compare beta
  ## values between localities
  samps <- ptrans(samps, grep('beta', colnames(samps)), 'pbeta')

  samps
}



merge_ensembles <- function(filenames, N)
{
  allsamps <- do.call(rbind, lapply(filenames, read1ensemble))

  ## Meh.  We can't use the thinning in metrosamp's getsamples function because
  ## we have to make some modifications to the samples.  Recreate the process
  ## here.
  ntot <- nrow(allsamps)
  i1 <- floor(ntot/N)
  i2 <- ntot - ntot%%N
  keep <- round(seq(i1, i2, length.out=N))

  allsamps[keep, , drop=FALSE]
}


get_betadist <- function()
{
  counties <- sort(unique(va_weekly_ntest_county$locality))

  fdir <- here::here('analysis', 'bayesian-filter-tests')

  betadists <- lapply(counties,
                      function(locality) {
                        fname <- file.path(fdir,
                                           paste0('filter-ensemble-candidates-',
                                                  locality, '.rds'))
                        if(!file.exists(fname)) {
                          return(NULL)
                        }
                        m <- getsamples(readRDS(fname))
                        betacol <- grep('beta', colnames(m))
                        m[ , betacol]
                      })
  names(betadists) <- counties

  betadists
}

save_harmonized_ensemble <- function(grand_ensemble, betadists)
{
  for(locality in names(betadists)) {
    if(!is.null(betadists[[locality]])) {
      ens <- grand_ensemble
      betacol <- which(colnames(ens) == 'pbeta')
      pbeta <- ens[ , betacol]
      betavals <- quantile(ecdf(betadists[[locality]]), pbeta)
      ens[ , betacol] <- betavals
      colnames(ens)[betacol] <- paste0('beta_', locality)

      filename <- paste0('filter-ensemble-',locality,'.rds')
      filename <- here::here('analysis','harmonized-ensemble', filename)
      saveRDS(ens, filename)
    }
  }
}
