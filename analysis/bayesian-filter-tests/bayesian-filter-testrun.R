parms <-
  c(beta = 0.33, A0 = 4, D0 = 4, Ts = 4, b0=10, b1=0.25, mask_effect = 0, import_cases=0)
vars <-
  c(S = 10000, E = 30, I = 0, Is = 0, R = 0)
tv0 <- seq(1,88)
obsdata <- CovMitigation:::get_obsdata()[[1]] %>% filter(locality=='FairfaxCounty')

run0 <- run_parmset(parms, vars, tv0)

istate <- c(run0[nrow(run0), ], parms)
istates <- do.call(rbind, rep(list(istate),100))

rslt1 <- bayesian_filter(istates, 'FairfaxCounty', 88, 95, obsdata)
rslt2 <- bayesian_filter(rslt1, 'FairfaxCounty', 95, 102, obsdata)
rslt3 <- bayesian_filter(rslt2, 'FairfaxCounty', 102, 109, obsdata)
rslt4 <- bayesian_filter(rslt3, 'FairfaxCounty', 109, 116, obsdata)
rslt5 <- bayesian_filter(rslt4, 'FairfaxCounty', 116, 123, obsdata)
rslt6 <- bayesian_filter(rslt5, 'FairfaxCounty', 123, 130, obsdata)
rslt7 <- bayesian_filter(rslt6, 'FairfaxCounty', 130, 137, obsdata)
rslt8 <- bayesian_filter(rslt7, 'FairfaxCounty', 137, 144, obsdata)
rslt9 <- bayesian_filter(rslt8, 'FairfaxCounty', 144, 151, obsdata)
rslt10 <- bayesian_filter(rslt9, 'FairfaxCounty', 151, 158, obsdata)

addprevalence <- function(rslt, pop) {
  rslt <- as.data.frame(rslt)
  rslt$fi <- (rslt$I + rslt$Is) / pop
  rslt
}

## We're not really running Fairfax County; we started with a population of 
## just 10k susceptible and 30 exposed
ffxpop <- 10030

rslts_by_time <- lapply(list(rslt1, rslt2, rslt3, rslt4, rslt5,
                             rslt6, rslt7, rslt8, rslt9, rslt10),
                        addprevalence, pop=ffxpop)


timevals <- seq(95, 158, 7)

statcols <- function(tbl, colname)
{
  probs=c(0.025, 0.5, 0.975)
  stats <- quantile(tbl[[colname]], probs, names=FALSE)
  names(stats) <- paste0(colname, c('lo','','hi'))
  stats
}

collate_results <- function(rsltlist, timevals)
{
  rlist <- lapply(seq_along(rsltlist),
                  function(i) {
                    t <- timevals[i]
                    rslt <- rsltlist[[i]]
                    probs <- c(0.025, 0.5, 0.975)
                    betastats <- statcols(rslt, 'beta')
                    impstats <- statcols(rslt, 'import_cases')
                    fistats <- statcols(rslt, 'fi')
                    row <- c(time=t, betastats, impstats, fistats)
                    m <- matrix(row, nrow=1)
                    d <- as.data.frame(m)
                    colnames(d) <- names(row)
                  })
}
final_rslt <- collate_results(rslts_by_time, timevals) %>% bind_rows()
