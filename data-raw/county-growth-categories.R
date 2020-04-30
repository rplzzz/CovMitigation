suppressPackageStartupMessages(library(dplyr))
library(tidyr)
library(readr)
library(NYTimesCOVID19)
library(here)

## Get the growth rates for the month of April
gr <- function(date, cases) {
  imax <- length(date)
  stopifnot(imax == length(cases))
  (log(cases[imax]) - log(cases[1])) / as.numeric(date[imax]-date[1])
}
aprildata <- 
  filter(cov19county, state == 'Virginia',
         !is.na(fips),
         date >= as.Date('2020-04-01'),
         date <= as.Date('2020-04-30')
         )
growthrates <- 
  group_by(aprildata, county, fips) %>%
  summarise(rate=gr(date, cases), cases=max(cases)) %>%
  ungroup() %>% 
  mutate(td=log(2)/rate, fips=as.integer(fips)) %>%
  arrange(fips)

## Translate the names to the ones we use internally in the model
fipscodes <- 
  read_csv(here('data-raw','VA_county_FIPScodes.csv'), col_types = 'ic') %>%
  rename(fips='FIPS_code', locality='Locality')

growthrates <- left_join(growthrates, fipscodes, by='fips')

## separate into categories
uhi_counties <- filter(growthrates, td <= 6.1) %>% select(fips, locality)
hi_counties <- filter(growthrates, td > 6.1, td <= 10) %>% select(fips, locality)
lo_counties <- filter(growthrates, td > 10, td <= 20) %>% select(fips, locality)
ulo_counties <- filter(growthrates, td > 20) %>% select(fips, locality)

### A few places don't yet have any observations, so put them in the same categories
### as nearby places.
## Martinsville -> Henry
## Bland -> Giles/Pulaski/Smyth
## Dickenson -> Buchanan & Wise are ultra-high, but I'm reluctant to throw a county
##              into that category without evidence
addtohi <- c('Martinsvillecity','BlandCounty', 'DickensonCounty')
## Bath -> Rockbridge/Alleghany
## Grayson -> Carroll/Wythe/GalaxCity  Smyth is also just under the cutoff for high.
addtolo <- c('BathCounty', 'GraysonCounty')               

hi_counties <- bind_rows(hi_counties, filter(fipscodes, locality %in% addtohi))
lo_counties <- bind_rows(lo_counties, filter(fipscodes, locality %in% addtolo))

usethis::use_data(uhi_counties, hi_counties, lo_counties, ulo_counties)
