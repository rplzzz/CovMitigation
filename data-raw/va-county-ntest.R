#### Parse the record of COVID-19 tests by district and assign to localities

library('here')
library('dplyr')
library('tidyr')

cumul_tests_by_district <- 
  readr::read_csv(here('data-raw','testing-by-district.csv'), 
                  col_types =  'ciiiiiiiiiiiiii')
cumul_tests_by_district$date <- lubridate::mdy(cumul_tests_by_district$date)

## Note we start this time series at min(cumul_tests_by_district$date) + 1 because
## we won't have the number of tests on that first date.
tottest <- 
  select(vacovdata::vacovdata, date, ntest = totalTestResultsIncrease,
         cumntest = totalTestResults) %>%
  arrange(date)
cumul_tests_by_district <- filter(cumul_tests_by_district,
                                  date >= min(tottest_meas$date)-1, date <= max(tottest_meas$date))
tottest_meas <- filter(tottest, 
                       date > min(cumul_tests_by_district$date), 
                       date <= max(cumul_tests_by_district$date))

ndate <- nrow(cumul_tests_by_district)  

ctbd <- as.matrix(dplyr::select(cumul_tests_by_district, -date))
tbd <- apply(ctbd, 2, diff)
stopifnot(nrow(tbd) == nrow(tottest_meas))

## There is one day (so far) on which the cumulative total inexplicably declines
## Assign this and any other such cases a single test (in fact there was a single
## case reported in the affected district that day)
tbd[tbd<0] <- 1

meas_district_total <- apply(tbd, 1, sum)
date_meas <- cumul_tests_by_district$date[seq(2,ndate)]

## subtract total measured from state total to get Rest of Virginia.
ROVA <- tottest_meas$ntest - meas_district_total

tests_by_district_meas <- cbind(date_meas, as.data.frame(tbd), ROVA) %>% rename(date=date_meas)

### Now deal with the dates prior to when they started reporting tests by district
tottest_unmeas <- filter(tottest, 
                         date <= min(cumul_tests_by_district$date),
                         !is.na(ntest))

dist_ratios <- ctbd[1, ,drop=FALSE] / as.numeric(tottest_unmeas[nrow(tottest_unmeas), 'cumntest'])
ROVA <- 1-sum(dist_ratios)
dist_ratios <- cbind(dist_ratios, ROVA)
ntest_unmeas <- matrix(tottest_unmeas$ntest, ncol=1) %*% dist_ratios
date <- tottest_unmeas$date
tests_by_district_unmeas <- cbind(date, as.data.frame(ntest_unmeas))

### And now VDH has stopped publishing this information.  I'm not really sure what
### the best approximation to make here is.  In theory, for dates before and after
### they started publishing the tests by district we should marginalize over the 
### unknown test allocation, but realistically we're not going to make that happen
### in a reasonable time frame.  Instead, we just take the cumulative totals as of
### the last data we have and keep those ratios going forward.
tottest_post_meas <- filter(tottest,
                            date > max(cumul_tests_by_district$date),
                            !is.na(ntest))
dist_ratios_post <- ctbd[nrow(ctbd), , drop=FALSE] / as.numeric(tottest_post_meas[1, 'cumntest'])
ROVA <- 1 - sum(dist_ratios_post)
dist_ratios_post <- cbind(dist_ratios_post, ROVA)
ntest_post <- matrix(tottest_post_meas$ntest, ncol=1) %*% dist_ratios_post
date <- tottest_post_meas$date
tests_by_district_post_meas <- cbind(date, as.data.frame(ntest_post))

tests_by_district <- rbind(tests_by_district_unmeas, tests_by_district_meas,
                           tests_by_district_post_meas)

### Splendid.  Now we need to allocate the tests to individual jurisdictions within 
### each district.  Start by loading the table of counties and districts.  It's
### pretty awesome that they give us the area of each county, but we won't need it.
vahd <- 
  readr::read_csv(here('data-raw','va-health-districts.csv'), 
                  col_types = 'iccnic') %>%
  select(-areaSqMile)  

## This table still has Bedford City and Clifton Forge City, which have since 
## merged into Bedford County and Alleghany county respectively.
ibedcty <- which(vahd$fips == 51515)
ibedcoun <- which(vahd$fips == 51019)
icliffrg <- which(vahd$fips == 51560)
ialleg <- which(vahd$fips == 51005)

vahd$population[ibedcoun] <- vahd$population[ibedcoun] + vahd$population[ibedcty]
vahd$population[ialleg] <- vahd$population[ialleg] + vahd$population[icliffrg]
vahd <- vahd[-c(ibedcty, icliffrg),]

## Now change all of the districts for which we don't have daily testing data
## to be part of the ROVA district.
vahd$district[!vahd$district %in% names(tests_by_district)] <- 'ROVA'
district_total_pop <- 
  group_by(vahd, district) %>%
  summarise(disttotpop = sum(population))
vahd <- 
  left_join(vahd, district_total_pop, by='district') %>% 
  mutate(distpopfrac=population/disttotpop) %>%
  select(fips, locality, district, population, distpopfrac)

## Now we have a table of tests by date and district, and a table of localities,
## districts, and population fractions.  A full join gives us all the information
## we need.
va_county_ntest <- 
  pivot_longer(tests_by_district, -date, names_to='district', values_to = 'district_ntest') %>%
  full_join(vahd, by='district') %>%
  mutate(ntest = district_ntest * distpopfrac) %>%
  select(date, fips, locality, district, population, ntest)

usethis::use_data(va_county_ntest, overwrite=TRUE)
