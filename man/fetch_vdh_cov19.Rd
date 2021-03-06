% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vdh.R
\name{fetch_vdh_cov19}
\alias{fetch_vdh_cov19}
\title{Download latest COVID-19 data from Virginia Department of Health}
\usage{
fetch_vdh_cov19()
}
\value{
Data frame with 5 columns:
\describe{
  \item{Report.Date}{Date cases were reported}
  \item{FIPS}{FIPS code for the locality the cases were reported in}
  \item{Locality}{Name of the locality the cases were reported in}
  \item{Health.District}{Name of the health district the locality belongs to}
  \item{Total.Cases}{Total number of cases (not clear whether this is for this day only, or 
  cumulative)}
}
}
\description{
Since this data changes frequently, we haven't included it as package data; instead,
this function will fetch the latest copy of the dataset from the VDH servers.
}
