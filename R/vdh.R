
#' Download latest COVID-19 data from Virginia Department of Health
#' 
#' Since this data changes frequently, we haven't included it as package data; instead,
#' this function will fetch the latest copy of the dataset from the VDH servers.
#' 
#' @return Data frame with 5 columns:
#' \describe{
#'   \item{Report.Date}{Date cases were reported}
#'   \item{FIPS}{FIPS code for the locality the cases were reported in}
#'   \item{Locality}{Name of the locality the cases were reported in}
#'   \item{Health.District}{Name of the health district the locality belongs to}
#'   \item{Total.Cases}{Total number of cases (not clear whether this is for this day only, or 
#'   cumulative)}
#' }
#' @export
fetch_vdh_cov19 <- function()
{
  vdh_url <- 'http://www.vdh.virginia.gov/content/uploads/sites/182/2020/03/VDH-COVID-19-PublicUseDataset-Cases.csv'
  vdh_data <- 
    utils::read.csv(vdh_url, stringsAsFactors = FALSE)
  ## First column name sometimes has odd characters in it.
  names(vdh_data)[1] <- 'Report.Date'
  stopifnot(setequal(names(vdh_data), 
                     c('Report.Date', 'FIPS', 'Locality', 'Health.District', 'Total.Cases')))
  vdh_data[['Report.Date']] <- lubridate::mdy(vdh_data[['Report.Date']])

  vdh_data
}
