library(CovMitigation)
library(deSolve)
library(tidyr)
library(ggplot2)

initvals <- c(S=999, I=1, R=0)
duration <- 7
beta0 <- as.numeric(2.5/(duration*initvals['S']))
params_const <- list(gamma=1/duration, beta=beta0)

rslt_const <- ode(initvals, seq(0,100), sir_equations, params_const)

params_step <- list(gamma=1/duration,
                 beta=data.frame(time=c(0,24), beta=c(beta0, beta0/2)))
rslt_step <- ode(initvals, seq(0, 100), sir_equations, params_step)

infections <- data.frame(time=rslt_const[,'time'], const_params=rslt_const[,'I'], 
                         step_params=rslt_step[,'I'])
infections <- pivot_longer(infections, -time, names_to = 'simulation type', values_to='Infections')

ggplot(data=infections, aes(x=time, y=Infections, color=`simulation type`)) +
  geom_line(size=1.25) +
  theme_bw(14)
