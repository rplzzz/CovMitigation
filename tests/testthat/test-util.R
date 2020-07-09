context('Utility functions')

library(dplyr)

test_that("Week aggregation works", {
  strtdate <- as.Date('2020-01-01')
  nweek <- 3
  firstday <- 4
  ndaywk <- 7
  x0 <- 10
  scl <- 7
  
  ## columns for the data
  time <- seq(0, firstday + ndaywk*nweek)
  date <- time + strtdate
  week <- ceiling((time-firstday)/ndaywk)
  x <- x0*exp(time/scl)
  y <- x
  
  df <- data.frame(date=date, time=time, week=week, x=x, y=y)
  
  weekending <- date[(time %% ndaywk) == (firstday %% ndaywk)]
  
  dfwkmean <- wkagg(weekending, df, c('time', 'x'))
  
  cmpdata_mean <- 
    as.data.frame(
      group_by(df, week)  %>%
        summarise(date=max(date), time=mean(time), x=mean(x), y=max(y)) %>%
        select(date, time, week, x, y)
    )
  
  expect_equal(dfwkmean, cmpdata_mean)
  
  ## Check that changing the aggregating function works
  dfwkmin <- wkagg(weekending, df, c('time','y'), aggfun=min)
  cmpdata_min <- 
    as.data.frame(
      group_by(df, week)  %>%
        summarise(date=max(date), time=min(time), x=max(x), y=min(y)) %>%
        select(date, time, week, x, y)
    )
  expect_equal(dfwkmin, cmpdata_min)
  
  ## Check that dates out of range are ignored.
  weekending <- c(as.Date('2019-01-01'), weekending)
  dfwkmean <- wkagg(weekending, df, c('time','x'))
  expect_equal(dfwkmean, cmpdata_mean)
  
  
})
