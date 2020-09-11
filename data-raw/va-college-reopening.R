#### Location and dates of colleges reopening in Virginia

library('readr')
library('tidyr')
library('dplyr')
library('here')

va_college_reopening <- read_csv(here('data-raw','va-college-reopening.csv'), col_types='ccccc') %>%
  rename(date=`Student Move-in date`) %>%
  mutate(date = lubridate::mdy(date)) %>%
  pivot_longer(c('locality1','locality2','locality3'), values_to = 'locality') %>%
  select(School, date, locality) %>%
  filter(!is.na(locality))

stopifnot(all(va_college_reopening$locality %in% vdhcovid::valocalities$locality))

usethis::use_data(va_college_reopening, overwrite=TRUE)
