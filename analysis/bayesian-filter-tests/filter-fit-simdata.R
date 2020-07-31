#### Fit a model to the base model.  We will start the filter fits at t=70

### 0. load the simulated datasets
source(here::here('analysis','bayesian-filter-tests','simulated-obs.R'))

### 1. Generate ensemble parameters and initial states.  
## parameters for the ensemble
Nens <- 100  # Ensemble size
R0 <- rnorm(Nens, 1.32, 0.1)
A0 <- rnorm(Nens, 4, 0.5)  # Should probably do these as gamma dist. to ensure nonnegative.
D0 <- rnorm(Nens, 4.5, 0.5)
Ts <- rnorm(Nens, 4, 1)
b_nt1 <- rlnorm(Nens, 4, 0.2)    # b value when nsamples == 1
b_nt100 <- rlnorm(Nens, 3.5, 0.2)   # b value when nsamples == 100
mask_effect <- rep(0, Nens)
import_cases <- rep(0, Nens)
## Convert these to the parameters used in our model
beta <- R0 / D0
b0 <- b_nt1
b1 <- (b0 - b_nt100)/log(100)
## create the parameter matrix
parms <- matrix(c(beta, A0, D0, Ts, b0, b1, mask_effect, import_cases), 
                nrow=Nens)
colnames(parms) <- c('beta','A0','D0','Ts','b0','b1','mask_effect','import_cases')
### 2. Run the parameter sets up to the start of the filtering.
## We don't really know when the infection started or how much initial exposure
## there was.  The dataset has 1/20 positive tests at t=21, so start with that
## as a guess, adjusting for our enrichment factor
obsdata <- obsdata_base
i0 <- which.max(obsdata$npos > 0)
t0 <- obsdata$time[i0]
fib <- obsdata$npos[i0] / obsdata$ntest[i0]
fi_odds <- fib/(1-fib) / (b0 - b1*log(obsdata$ntest[i0]))
fi <- fi_odds / (1+fi_odds)
E0 <- fi * obsdata_base$population
S0 <- obsdata_base$population - E0
vars0 <- cbind(t=t0, S=S0, E=E0, I=0, Is=0, R=0)

t1 <- 70
timevals <- seq(t0,t1)
vars1 <- t(sapply(seq(1,Nens), 
                   function(i) {
                     x <- run_parmset(parms[i,], vars0[i,], timevals)
                     x[nrow(x),]
                   }))
state1 <- cbind(vars1, parms)
fit_base <- fit_filter(state1, obsdata, t1, max(obsdata$time), history = TRUE)

