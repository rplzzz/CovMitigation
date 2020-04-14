#### Documentation for datasets

#' Census data for Virginia localities
#' 
#' Dataset with total population and age breakdown for Virginia localities (counties,
#' county-equivalents, and the entire commonwealth).
#' 
#' There were no entries for Suffolk City, Virginia Beach City, or Williamsburg City in the dataset, 
#' so we looked up the total
#' population and estimated the age and sex breakdowns by applying the category
#' ratios for Norfolk, a nearby city that did have data.
#' 
#' Winchester City also did not have data, so we used the ratios for surrounding 
#' Frederick County.
#' 
#' Obviously
#' it would be better to get the real data.
#' 
#' @format Data frame with 8 columns
#' \describe{
#' \item{Locality}{Name of the locality}
#' \item{Total}{Total population}
#' \item{Age0to5}{Population ages 0--5}
#' \item{Age5to19}{Population ages 5--19}
#' \item{Age20to39}{Population ages 20--39}
#' \item{Age50to59}{Population ages 50--59}
#' \item{Age60to74}{Population ages 60--74}
#' \item{Age75plus}{Population over age 75}
#' }
#' 
#' @source Weldon Cooper Center (verify?)
"vaMyAgeBands"

#' Market fractions for Virginia localities
#' 
#' Dataset giving the market fractions by locality (county or county-equivalent) for UVA and VCU, 
#' for all DRGs and for the Respiratory DRG
#' 
#' @format Data frame with 3 columns
#' \describe{
#' \item{Locality}{Name of the locality}
#' \item{TypeAMCDRG}{Synthetic key combining the type (always 'fraction'), DRG (either 'All' or 'Resp'),
#' and medical center (either 'UVA' or 'VCU')}
#' \item{marketShare}{The medical center's market share in that locality, for that DRG, expressed as a 
#' fraction}
#' }
#' @source (not given)
"marketFractionFinal"

#' Market fraction and FIPS code for Virginia localities with high UVA market fractions
#' 
#' Same as \code{\link{marketFractionFinal}}, except includes an additional column for FIPS
#' code and hs been filtered to just localities with UVA market shares >= 20%
"sampleCounties"

#' Virginia Department of Health COVID-19 surveillance data
#' 
#' Dataset with reported COVID-19 cases by date and locality.  It is not yet clear whether the
#' "Total Cases" reported are cumulative or new on that day.
#' 
#' @format Data frame with 5 columns:
#' \describe{
#' \item{Report.Date}{Date of the report}
#' \item{FIPS}{FIPS code for the locality}
#' \item{Locality}{Name of the locality}
#' \item{Health.District}{Health district the locality belongs to}
#' \item{Total.Cases}{Number of COVID-19 cases reported.  It's not documented in
#' the original source whether this is cumulative or for that day only, but the data
#' matches the figures given on the agency's website, which appears to be cumulative.}
#' }
#' @source Virginia Department of Health:\cr
#' \code{http://www.vdh.virginia.gov/content/uploads/sites/182/2020/03/VDH-COVID-19-PublicUseDataset-Cases.csv}
"vdh_covid19"

#' Reprocessed data from VA Department of Health
#' 
#' This dataset contains the VDH COVID-19 data manually collected by the UVAHS 
#' data science team.  It contains some data not present in the VDH daily download,
#' such as the number of tests performed; however, data is only presented at the 
#' commonwealth and TJ health district aggregate levels.
#' 
#' @format Data frame with 10 columns:
#' \describe{
#' \item{date}{Report date}
#' \item{vaNewCases}{New cases reported across the commonwealth on that date.}
#' \item{vaCumCases}{Cumuative cases reported across the commonwealth.  This is 
#' every instance of infection ever reported, irrespective of whether the patient is 
#' still infected.}
#' \item{tjNewCases}{New cases reported in the Thomas Jefferson health district.}
#' \item{tjCumCases}{Cumulative cases reported in the Thomas Jefferson health district.}
#' \item{nhosp}{Number of COVID-19 cases admitted to hospitals (commonwealth total).}
#' \item{ntest}{Number of COVID-19 tests performed (commonwealth total).}
#' \item{ftest}{Fraction of total population tested on this date.}
#' \item{ntest_cum}{Cumulative number of tests performed.}
#' \item{posTestFrac}{Fraction of tests with positive results (=vaNewCases/ntest).}
#' \item{vapop}{Total population of Virginia}
#' }
#' @source Virginia Department of Health:  \code{http://www.vdh.virginia.gov/coronavirus/}
"uvads_covid19"

#' Date of first report for VA counties
#' 
#' Date of the first reported COVID-19 case for each Virginia county.  Counties
#' that have not yet had a reported case have a date of \code{NA}.
#' 
#' @format Data frame with 4 columns
#' \describe{
#' \item{FIPS}{FIPS code for the county.}
#' \item{firstReport}{Date of the first reported COVID-19 case in the county.}
#' \item{firstDay}{Number of days since 1Jan2020 for the firstReport}
#' }
#' @source New York Times COVID-19 dataset \code{https://github.com/nytimes/covid-19-data}
"va_county_first_case"

#' Counts of various categories of COVID cases
#' 
#' @format Data frame with 10 columns
#' \describe{
#' \item{date}{Date of observation}
#' \item{COVID}{Total number of COVID-19 cases}
#' \item{Admits}{Number of COVID-19 cases admitted that day}
#' \item{DCs}{Discharges}
#' \item{ICU}{Number of COVID-19 cases in ICU}
#' \item{ICU in}{Number of new ICU cases}
#' \item{ICU out}{Number of cases leaving the ICU}
#' \item{Vent}{Number of cases on IMV}
#' \item{vent in}{Number of new IMV cases}
#' \item{vent out}{Number of IMV cases leaving}
#' }
"uva_covid_count"
