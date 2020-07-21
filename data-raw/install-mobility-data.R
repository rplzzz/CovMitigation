library(readr)
library(dplyr)
library(CovMitigation)

dataurl <- 'https://www.gstatic.com/covid19/mobility/Global_Mobility_Report.csv'

colspec <- 'ccccciDnnnnnn'

rawdata <- readr::read_csv(dataurl, col_types=colspec) %>%
  filter(country_region_code=='US',
         sub_region_1=='Virginia',
         !is.na(census_fips_code)) %>%
  select(date, fips=census_fips_code, county=sub_region_2, 
         retail = retail_and_recreation_percent_change_from_baseline,
         grocery = grocery_and_pharmacy_percent_change_from_baseline,
         parks = parks_percent_change_from_baseline,
         transit = transit_stations_percent_change_from_baseline,
         work = workplaces_percent_change_from_baseline,
         home = residential_percent_change_from_baseline
  )

## Add the version of the locality identifiers used in state data
locality_names <- select(vdhcovid::valocalities, fips, locality)
va_mobility_temp <- left_join(rawdata, locality_names, by='fips')

## Fill in missing data using loess smoothing
va_mobility_daily <- 
  group_by(va_mobility_temp, fips, county, locality) %>%
  group_map(
    function(df, group) {
      df <- tidyr::complete(df, date=seq(min(df$date), max(df$date), by=1))
      t <- as.numeric(df[['date']] - df[['date']][1])
      rslt <- as.list(group)
      rslt[['date']] <- df[['date']]
      for(col in c('retail','grocery','parks','transit','work','home')) {
        y <- df[[col]]/100      ## Convert percentage to fraction.
        miss <- is.na(y)
        if(sum(!miss) < 20) {
          rslt[[col]] <- NA_real_
        }
        else {
          lo <- loess(y ~ t, data=NULL, na.action = na.omit, span=0.25 )
          y[miss] <- predict(lo, t[miss])
          rslt[[col]] <- y
        }
      }
      tibble::as_tibble(rslt)
    }
  ) %>%
  bind_rows() %>%
  mutate(t = as.numeric(date - as.Date('2020-01-01')),
         home = -home           # Give home index the same sign convention as other indices
  )

daily2weekly <- function(dailydata) {
  strt <- as.Date('2019-12-30')
  t <- as.numeric(dailydata$date - strt)
  dailydata$week <- as.integer(floor(t/7))
  weeklydata <- 
    group_by(dailydata, fips, locality, week) %>%
    summarise(date = week[1]*7 + strt, retail=mean(retail), grocery=mean(grocery), parks=mean(parks), work=mean(work),
              home=mean(home)) %>%
    ungroup() %>%
    mutate(t = as.numeric(date - as.Date('2020-01-01')))
  
  sentinels <- lapply(unique(weeklydata$locality), 
                      function(loc) {
                        ld <- filter(weeklydata, locality==loc)
                        sentinel <- ld[nrow(ld),]
                        for(col in c('retail', 'grocery', 'parks', 'work','home')) {
                          if(all(is.na(ld[[col]]))) {
                            sentinel[1,col] <- 0
                          }
                          else {
                            good <- ld[[col]][!is.na(ld[[col]])]
                            sentinel[1,col] <- good[length(good)]
                          }
                        }
                        sentinel[1,'t'] <- 1e6
                        sentinel[1,'date'] <- 1e6 + as.Date('2020-01-01')
                        sentinel
                      }) %>%
    bind_rows()
  
  bind_rows(sentinels, weeklydata) %>%
    arrange(locality, t)
}

va_mobility_weekly <- daily2weekly(va_mobility_daily)


### Create future mobility scenarios
## quadratic return to normal mobility on 2020-09-01
future_start <- max(va_mobility_daily$date)
future_end <- as.Date('2020-09-01')
jan01 <- as.Date('2020-01-01')
fdate <- seq(future_start, future_end, by=1)
ftime <- as.numeric(fdate - jan01)
x1 <- min(ftime)
x2 <- max(ftime)
va_future_mobility_daily <- 
  lapply(unique(va_mobility_daily$locality),
         function(loc) {
           localmob <- dplyr::filter(va_mobility_daily, locality==loc)
           imax <- which.max(localmob$date)
           datacols <- c('retail', 'grocery', 'parks', 'transit', 'work', 'home')
           fmd <- tibble::tibble(fips=localmob$fips[1], county=localmob$county[1],
                                 locality=loc, date=fdate)
           for(col in datacols) {
             gooddata <- localmob[!is.na(localmob[[col]]), ]
             if(nrow(gooddata) == 0) {
               fmd[[col]] <- 0
             }
             else {
               imax <- which.max(gooddata$date)
               y1 <- gooddata[[col]][imax]
               y2 <- 0.25
               if(y1 >= y2) {
                 fmd[[col]] <- y1
               }
               else {
                 b <- (y1-y2) / (-x1^2/(2*x2) + x1 - x2/2)
                 c <- -b*x2/2 + y2
                 a <- -b/(2*x2)
                 fmd[[col]] <- a*ftime^2 + b*ftime + c
               }
             }
           }
           
           dplyr::bind_rows(localmob, fmd)
         }) %>%
  dplyr::bind_rows()

va_future_mobility_weekly <- daily2weekly(va_future_mobility_daily)

usethis::use_data(va_mobility_daily, va_mobility_weekly, 
                  va_future_mobility_daily, va_future_mobility_weekly,
                  overwrite=TRUE)
