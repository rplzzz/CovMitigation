## Update the Virginia Department of Health COVID-19 dataset

## This is the version of the data that the UVAHS data science team has been
## collecting manually from the VDH website
library(readr)
library(dplyr)
library(here)


vdh_new <- read_csv(here('data-raw', 'Virginia_Dept_Health_COVID_Data.csv'),
                    col_types = 'ciiiiinnniii')

uvads_covid19 <-
  filter(vdh_new, !is.na(Date)) %>%
  rename(posTestFrac=daily_pos_pct_of_tests, date=Date, vaNewCases=VA_New_Confirmed_Cases,
         vaCumCases=VA_Cumulative_Confirmed_Cases, tjNewCases=TJ_Confirmed_Cases,
         tjCumCases=TJ_Cumulative_Cases, nhosp=num_hospitalizations, ntest_cum=cum_num_tests,
         ntest=daily_num_tests,
         vapop=va_population)


## Add a new tests column
#uvads_covid19$ntest <- c(NA, diff(uvads_covid19$ntest_cum))
uvads_covid19 <-
  mutate(uvads_covid19, date=lubridate::mdy(date), ftest=ntest/vapop) %>%
  select(date, vaNewCases, vaCumCases, tjNewCases, tjCumCases, nhosp, 
         ntest, ntest_cum, ftest, posTestFrac, vapop)


usethis::use_data(uvads_covid19, overwrite=TRUE)


## Hospitalization dataset
uva_covid_count <- readr::read_csv(here('data-raw','COVID_counts_uva.csv'), 
                                   col_types = 'ciiiiiiiii')
uva_covid_count$date <- lubridate::mdy(uva_covid_count$date)
usethis::use_data(uva_covid_count, overwrite=TRUE)
