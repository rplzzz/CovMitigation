library(CovMitigation)

set.seed(867-5309)
#### Set up baseline scenario for testing
base_parms <-
  c(beta = 0.25, A0 = 4, D0 = 4.5, Ts = 4, b0=50, b1=4, mask_effect = 0, import_cases=0)
poptot <- 100000
E0 <- 200
vars <- c(t=0, S = poptot-E0, E = E0, I = 0, Is = 0, R = 0)
tpast <- seq(1, 175)

base_run <- run_parmset(base_parms, vars, tpast)

obsdata_base <- simobs(base_run, base_parms[c('b0','b1')])

#### Ramp up scenario:  beta increases to 0.4 between t=60 and t=120
t1 <- 154; t2 <- 174
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

obsdata_ramp <- simobs(ramp_run, ramp_parms[c('b0','b1')])
