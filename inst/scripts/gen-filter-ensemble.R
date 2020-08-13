library(CovMitigation)
library(vdhcovid)
library(dplyr)
library(foreach)
library(metrosamp)

counties <- sort(unique(vdhcovid::va_weekly_ntest_county$locality))
high_beta_counties <-
  c('FairfaxCounty', 'Fairfaxcity','FallsChurchcity', 'ArlingtonCounty',
    'Alexandriacity',
    'PrinceWilliamCounty', 'Manassascity', 'ManassasParkcity',
    'Richmondcity', 'HenricoCounty', 'ChesterfieldCounty')

gen_ensemble <- function(icounty)
{
  foreach::registerDoSEQ()

  county <- counties[icounty]

  common_parms <- c(D0=5.5, A0=2, I0=30, Ts=5.5, mask_effect=-0.5, b0=100, b1=1)
  beta <-
    if(county %in% high_beta_counties) {
      0.25
    } else {
      0.2
    }
  p0 <- c(beta, common_parms)
  names(p0) <- c(paste0('beta_',county), names(common_parms))

  lpost <- gen_post(fixed_parms=p0, maxdate = as.Date('2020-06-30'))

  message('county = ', county)

  pstrt <- p0
  ctrl <- list(fnscale=-1, maxit=2000)
  opt_uncons <- optim(pstrt, lpost, control=ctrl)

  set.seed(867-5309)

  scl <- c(0.01, 0.05, 0.5, 0.05, 0.05, 0.05, 1.0, 0.05)
  names(scl) <- names(p0)
  ms <- metrosamp(lpost, opt_uncons$par, 50000, 1, scl)

  outfilename <- paste0('filter-ensemble-candidates-',county, '.rds')
  saveRDS(ms, outfilename)
  message('Output file is ',outfilename)
  message('FIN.')
  outfilename
}


