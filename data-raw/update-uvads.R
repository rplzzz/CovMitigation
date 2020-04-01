## Update the Virginia Department of Health COVID-19 dataset

## This is the version of the data that the UVAHS data science team has been
## collecting manually from the VDH website
library(readr)
library(dplyr)
library(lubridate)
library(here)


vdh_new <- read_csv(here::here('data-raw', 'Virginia_Dept_Health_COVID_Data.csv'),
                    col_types = 'ciiiiinniinnn')
uvads_covid19 <-
  filter(vdh_new, !is.na(Date)) %>%
  rename(posTestFrac=pos_pct_of_tests, date=Date, vaNewCases=VA_New_Confirmed_Cases,
         vaCumCases=VA_Cumulative_Confirmed_Cases, tjNewCases=TJ_Confirmed_Cases,
         tjCumCases=TJ_Cumulative_Cases, nhosp=num_hospitalizations, ntest=num_tests,
         vapop=va_population) %>%
  mutate(date=mdy(date), ftest=ntest/vapop) %>%
  select(date, vaNewCases, vaCumCases, tjNewCases, tjCumCases, nhosp, 
         ntest, ftest, posTestFrac, vapop)


usethis::use_data(uvads_covid19, overwrite=TRUE)
