#### Fit a model to the base model.  We will start the filter fits at t=70

### 0. load the simulated datasets
source(here::here('analysis','bayesian-filter-tests','simulated-obs.R'))

### 1. Generate ensemble parameters and initial states.  
## parameters for the ensemble:  Start with an initial guess and optimize
pinit <- c(beta=0.25, A0=3, D0=4, Ts=4, b0=50, b1=1, import_cases=0, I0=30)
lpost_base <- gen_simobs_posterior(obsdata_base)
opt_filter_base <- optim(pinit, lpost_base, control=list(fnscale=-1, maxit=2000))

## data for diagnostic plots
popt_base <- opt_filter_base$par
popt_base['mask_effect'] <- 0
vars <- c(time=0, S=obsdata_base$population[1]-popt_base[['I0']], E=popt_base[['I0']],
          I=0, Is=0, R=0)
opt_filter_base_df <- as.data.frame(run_parmset(popt_base, vars, obsdata_base$time))
opt_filter_base_df$fi <- (opt_filter_base_df$I + opt_filter_base_df$Is) / obsdata_base$population[1]

Nens <- 100  # Ensemble size
Nsamp <- 10000 # Number of MCMC samples to use to generate the ensemble
scl <- c(beta=0.01, A0=0.25, D0=0.25, Ts=0.25, b0=5, b1=0.1, import_cases=1, I0=5)
ms_base <- metrosamp::metrosamp(lpost_base, opt_filter_base$par, Nsamp, 1, scl)
parms <- cbind(metrosamp::getsamples(ms_base, Nens), mask_effect=0)
## Rearrange the columns to the order we used before.
parms <- parms[ ,c('beta','A0','D0','Ts','b0','b1', 'I0', 'mask_effect','import_cases')]

### 2. Run the parameter sets up to the start of the filtering.
make_init_state <- function(parms, obsdata, tstrt, Nens) {
  vars1 <- t(sapply(seq(1,Nens), 
                   function(i) {
                     x <- initialize_parmset(parms[i,], obsdata, tstrt)
                     x[nrow(x),]
                   }))
  cbind(vars1, parms)
}

t1 <- 70
state_base <- make_init_state(parms, obsdata_base, t1, Nens)
fit_base <- fit_filter(state_base, obsdata_base, t1, max(obsdata_base$time), history = TRUE)
plots_base <- filter_model_diagnositc_plots(fit_base$modelfit, base_run)

## 2a. Run the filtering on the ramped beta data set
## TODO:  refactor this so that we can put in any data w want
state_ramp <- make_init_state(parms, obsdata_ramp, t1, Nens)
fit_ramp <- fit_filter(state_ramp, obsdata_ramp, t1, max(obsdata_ramp$time), history=TRUE)
plots_ramp <- filter_model_diagnositc_plots(fit_ramp$modelfit, ramp_run)

## 2b. Run on the import dataset (not run because the results appear to be mostly
## the same as the ramp run.)

## 2c. Run the filtering on the point-prevalence dataset
state_pp <- make_init_state(parms, obsdata_pp, t1, Nens)
fit_pp <- fit_filter(state_pp, obsdata_pp, t1, max(obsdata_pp$time), history=TRUE)
plots_pp <- filter_model_diagnositc_plots (fit_pp$modelfit, pp_run)
