library('CovMitigation')
library('vdhcovid')
library('doParallel')

set.seed(867-5309)
registerDoParallel(2)

int2locality <- function(i) 
{
  locs <- sort(unique(va_weekly_ntest_county$locality))
  locs[i]
}

run_filter_fit_loc <- function(i, N=100)
{
  ## Assume that the directory for filter ensemble candidates is a subdir of the
  ## directory we are running from
  inputdir <- normalizePath('./ensemble-candidates')
  i <- as.integer(i)
  locality <- int2locality(i)
  N <- as.integer(N)
  fit <- filter_fit_locality(locality, inputdir, N)
  
  outfilename <- paste0('filter-fit-', locality, '.rds')
  saveRDS(fit, outfilename)
  message('FIN.')
}