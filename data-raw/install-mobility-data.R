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

wkdates <- unique(vdhcovid::vaweeklytests[['date']])

va_mobility_weekly <- 
  group_by(va_mobility_daily, fips, locality) %>%
  group_map(function(df, group) {
    wkagg(wkdates, df, c('retail','grocery','parks','transit', 'work', 'home'))
  },
  keep=TRUE) %>%
  bind_rows()

usethis::use_data(va_mobility_daily, va_mobility_weekly, overwrite=TRUE)
