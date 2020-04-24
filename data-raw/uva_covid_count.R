library(here)

uva_covid_count <- readr::read_csv(here('data-raw','COVID_counts_uva.csv'),
col_types = 'ciiiiiiiii')
uva_covid_count$date <- lubridate::mdy(uva_covid_count$Day)
uva_covid_count$time <- as.numeric(uva_covid_count$date - as.Date('2020-01-01'))
usethis::use_data(uva_covid_count, overwrite=TRUE)

