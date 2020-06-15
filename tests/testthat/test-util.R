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
  
  cmpdata <- 
    as.data.frame(
      group_by(df, week)  %>%
        summarise(date=max(date), time=mean(time), x=mean(x), y=max(y)) %>%
        select(date, time, week, x, y)
    )
  
  expect_equal(dfwkmean, cmpdata)
  
  ## Check that changing the aggregating function works
  dfwkmean <- wkagg(weekending, df, c('time','y'), aggfun=min)
  cmpdata <- 
    as.data.frame(
      group_by(df, week)  %>%
        summarise(date=max(date), time=min(time), x=max(x), y=min(y)) %>%
        select(date, time, week, x, y)
    )
  expect_equal(dfwkmean, cmpdata)
  
  ## Check that we get an error if we ask for a date that isn't in the data frame
  weekending <- c(as.Date('2019-01-01'), weekending)
  expect_error(dfwkmean <- wkagg(weekending, df, c('time','y')))
  
  
})
