#### Documentation for datasets

#' Census data for Virginia localities
#' 
#' Dataset with total population and age breakdown for Virginia localities (counties,
#' county-equivalents, and the entire commonwealth).
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
