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
#' @source New York Times COVID-19 dataset \url{https://github.com/nytimes/covid-19-data}
"va_county_first_case"

#' Counts of various categories of COVID cases
#' 
#' @format Data frame with 10 columns
#' \describe{
#' \item{date}{Date of observation}
#' \item{time}{Days since 2020-01-01 for observation}
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


#' County growth categories
#' 
#' COVID-19 growth across the Commonwealth of Virginia has been very heterogeneous.
#' Over the month of April doubling times have ranged from under 4 days to nearly
#' 40 days.  The growth statistics over April have been used to sort the counties
#' into four categories.  In the model, each of these categories will get its own
#' growth rate.
#' 
#' @format Data frame with 2 columns
#' \describe{
#' \item{fips}{FIPS code}
#' \item{locality}{Name, using the vdh conventions}
#' }
#' @name growth_categories
NULL

#' @describeIn growth_categories High growth rate counties
"hi_counties"

#' @describeIn growth_categories Low growth rate counties
"lo_counties"

#' @describeIn growth_categories Ultra-high growth rate counties
"uhi_counties"

#' @describeIn growth_categories Ultra-low growth rate counties
"ulo_counties"

#' UVAHS Covid admissions by date and locality of origin
#' 
#' COVID cases at the HS tabulated by date of admission and locality.
#' 
#' @format Data frame with 4 columns
#' \describe{
#' \item{date} Date of admission
#' \item{fips} FIPS code for the locality of origin
#' \item{locality} Name of locality of origin
#' \item{hospitalizations} Number of hospitalizations
#' }
#' 
#' @source UVAHS records
"hospitalizations"
