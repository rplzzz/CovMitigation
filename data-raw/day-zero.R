## Dataset of date of first observation ov COVID-19 case in each VA locality
library(dplyr)
library(tidyr)
library(readr)
library(vdhcovid)
library(here)

va_county_first_case <-
  filter(vadailycases, cases > 0) %>%
  group_by(fips, locality) %>%
  summarise(firstReport = min(date)) %>%
  mutate(firstDay = as.numeric(firstReport - as.Date('2020-01-01'))) %>%
  rename(FIPS=fips, Locality=locality)

## For some reason, Bristol is in the list twice.  Fix the dupe.
va_county_first_case <- unique(va_county_first_case)

usethis::use_data(va_county_first_case, overwrite=TRUE)
