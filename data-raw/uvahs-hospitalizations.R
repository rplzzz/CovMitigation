#### Process data for UVA hosp admissions
library(here)
library(dplyr)
library(RODBC)

## Read queries from an SQL file and split into individual queries
## Assumes there is a tag "QSPLT" at the start of each query 
sqlread <- function(filename)
{
  sql <- readLines(filename)
  qstrt <- grep('^--QSPLT', sql)
  nq <- length(qstrt)
  qend <- c(qstrt[2:nq]-1, length(sql))
  
  sapply(seq_along(qstrt), function(i) {paste(sql[qstrt[i]:qend[i]], collapse='\n')})
}

queries <- sqlread(here('data-raw', 'covid-sql-queries', 'COVID_Positive_Hospitalizations.sql'))

DEFAULT_DB_CON <- 'driver={SQL Server};server=HSTSARTDM;database=DS_HSDM_COVID;trusted_connection=true'
dbcon <- odbcDriverConnect(DEFAULT_DB_CON)

## Assumption: all of the queries except the last one are selecting into temp
## tables on the server.  The last query result is the only one with data in it.
for(i in seq_along(queries)) {
  qrslt <- sqlQuery(dbcon, queries[i], stringsAsFactors = FALSE)
}
close(dbcon)

## We need to convert the Locality column in the data to the VDH formatting.  Easiest
## way to do this is to convert to FIPS
fipstbl <- readr::read_csv(here('data-raw', 'covid-sql-queries', 'Locality_FIPS.csv'), col_types = 'ci')
## All of the locality names in the result are in upper case
fipstbl$Locality <- toupper(fipstbl$Locality)
qrslt <- inner_join(qrslt, fipstbl, by='Locality') %>% rename(fips=FIPS)
## join in the vdh-style representation of locality
qrslt <- inner_join(qrslt, select(vdhcovid::valocalities, fips, locality), by='fips')

## We need to subtract off the time from symptom onset to hospitalization, so that
## the date reflects the approximate time the patient became symptomatic
tth <- 5
hospitalizations <-
  select(qrslt, date=Adm_Dtm, fips, locality) %>%
  mutate(date=as.Date(date)-5) %>%
  group_by(date, fips, locality) %>%
  summarise(hospitalizations = n()) %>%
  ungroup() %>%
  filter(date >= as.Date('2020-03-01')) %>%
  arrange(date)

## Aggregate daily values to weekly
strt <- as.Date('2019-12-30')              # Last Monday before 2020-01-01
wk0date <- as.Date('2020-01-05')           # Date for week==0.
t <- hospitalizations$date - strt
hospitalizations$week <- as.integer(floor(t/7))
hospitalizations <- 
  group_by(hospitalizations, week, fips, locality) %>%
  summarise(date=max(week)*7 + wk0date, hospitalizations = sum(hospitalizations)) %>%
  ungroup() %>%
  arrange(date)

#### Fill in missing date/locality combinations with zeros
## get all combinations of week and fips
weeks <- seq(min(hospitalizations$week), max(hospitalizations$week))
weekfips <- tidyr::expand_grid(week=weeks, fips=unique(vdhcovid::valocalities$fips))
wkdate <- unique(select(hospitalizations, week, date))
fipsloc <- unique(select(vdhcovid::valocalities, fips, locality))
## get the locality and date variables
weekfips <- 
  left_join(weekfips, wkdate, by='week') %>%
  left_join(fipsloc, by='fips')

## NA values are absent from the original dataset, so replace with zero
hospitalizations <-
  left_join(weekfips,hospitalizations, by=c('week','fips','date','locality')) %>%
  tidyr::replace_na(list(hospitalizations=0))

usethis::use_data(hospitalizations, overwrite=TRUE)
