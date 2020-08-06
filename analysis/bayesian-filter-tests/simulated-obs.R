library(CovMitigation)

set.seed(867-5309)
#### Set up baseline scenario for testing
base_parms <-
  c(beta = 0.25, A0 = 4, D0 = 4.5, Ts = 4, b0=50, b1=4, mask_effect = 0, import_cases=0)
poptot <- 100000
E0 <- 200
vars <- c(t=0, S = poptot-E0, E = E0, I = 0, Is = 0, R = 0)
tpast <- seq(1, 182)

base_run <- run_parmset(base_parms, vars, tpast)

obsdata_base <- simobs(base_run, base_parms[c('b0','b1')])

#### Ramp up scenario:  beta increases to 0.4 between t=60 and t=175
t1 <- 154; t2 <- 175
betabase <- 0.33; betamax <- 0.4
ramp_run <- base_run      # run is the same as baseline up through t=t1
ramptimes <- seq(t1, t2)
irows <- match(ramptimes, base_run[ , 'time'])
ramp_parms <- base_parms
for(i in seq_along(ramptimes)) {
  t <- ramptimes[i]
  irow <- irows[i]
  beta <- betabase + (t-t1+1)/(t2-t1+1) * (betamax-betabase)
  ramp_parms[['beta']] <- beta
  newrun <- run_parmset(ramp_parms, ramp_run[irow-1,], c(t-1, t))
  ramp_run[irow, ] <- newrun[1, ]
}
postramp_times <- seq(t2+1, max(tpast))
irowspost <- match(postramp_times, base_run[ , 'time'])
ramp_parms[['beta']] <- betamax
postrun <- run_parmset(ramp_parms, ramp_run[max(irows), ], c(t2, postramp_times))
ramp_run[irowspost, ] <- postrun


obsdata_ramp <- simobs(ramp_run, ramp_parms[c('b0','b1')])

#### Import surge scenario:  daily imports increase from 0 to 100 per day between
#### t=60 and t=175
t1 <- 154; t2 <- 175
importbase <- 0; importmax <- 100
import_run <- base_run      # run is the same as baseline up through t=t1
importtimes <- seq(t1, t2)
irows <- match(importtimes, base_run[ , 'time'])
import_parms <- base_parms
for(i in seq_along(importtimes)) {
  t <- importtimes[i]
  irow <- irows[i]
  imp <- importbase + (t-t1+1)/(t2-t1+1) * (importmax-importbase)
  import_parms[['import_cases']] <- imp
  newrun <- run_parmset(import_parms, import_run[irow-1,], c(t-1, t))
  import_run[irow, ] <- newrun[1, ]
}
postimport_times <- seq(t2+1, max(tpast))
irowspost <- match(postimport_times, base_run[ , 'time'])
import_parms[['import_cases']] <- importmax
postrun <- run_parmset(import_parms, import_run[max(irows), ], c(t2, postimport_times))
import_run[irowspost, ] <- postrun


obsdata_import <- simobs(import_run, import_parms[c('b0','b1')])


#### Point prevalence scenario:  Actual prevalence is as for base scenario, but
#### measured results are inflated by point prevalence results at t=140 and t=147
pp_run <- base_run
obsdata_pp <- obsdata_base
irows <- which(obsdata_pp$time %in% c(140,147))
obsdata_pp$ntest[irows] <- obsdata_pp$ntest[irows] + 100
obsdata_pp$npos[irows] <- obsdata_pp$npos[irows] + rbinom(2, 100, 0.95)

