## Dataset of date of first observation ov COVID-19 case in each VA locality
library(dplyr)
library(tidyr)
library(readr)
library(NYTimesCOVID19)
library(here)

vacounties <- read_csv(here::here('data-raw', 'VA_county_FIPScodes.csv'), col_types = 'ic') %>%
  rename(FIPS = FIPS_code) %>%
  drop_na()

county_first_case <-
  filter(cov19county, state=='Virginia', cases > 0) %>%
  rename(FIPS=fips) %>%
  mutate(FIPS=as.integer(FIPS)) %>%
  group_by(FIPS) %>%
  summarise(firstReport = min(date))

va_county_first_case <- left_join(vacounties, county_first_case, by='FIPS')
va_county_first_case$firstDay <- as.numeric(va_county_first_case$firstReport - as.Date('2020-01-01'))

## For some reason, Bristol is in the list twice.  Fix the dupe.
va_county_first_case <- unique(va_county_first_case)

usethis::use_data(va_county_first_case, overwrite=TRUE)
