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

hospitalizations <-
  select(qrslt, date=Adm_Dtm, fips, locality) %>%
  mutate(date=as.Date(date)) %>%
  group_by(date, fips, locality) %>%
  summarise(hospitalizations = n()) %>%
  ungroup() %>%
  arrange(date)

usethis::use_data(hospitalizations, overwrite=TRUE)
