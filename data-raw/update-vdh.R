## Update the Virginia Department of Health COVID-19 dataset

vdh_old <- CovMitigation::vdh_covid19
vdh_new <- CovMitigation::fetch_vdh_cov19()

vdh_nodup <- dplyr::anti_join(vdh_new, vdh_old, by=c('Report.Date', 'FIPS'))

if (nrow(vdh_nodup) == 0) {
  odatemax <- max(vdh_old$Report.Date)
  ndatemax <- max(vdh_new$Report.Date)
  message('No new data in current report.  Most recent date old:  ', odatemax,
          '  Most recent date new:  ', ndatemax)
} else {
  newdates <- paste(unique(vdh_nodup$Report.Date), collapse=', ')
  vdh_covid19 <- dplyr::bind_rows(vdh_old, vdh_nodup)
  ## We may not necessarily want to update the package data every time we get new data
  ## (especially if we are running on a cron job), so just save the results 
  #usethis::use_data(vdh_covid19, overwrite = TRUE)
  filename <- paste0('vdh-data-',as.Date(lubridate::now()), '.rds')
  saveRDS(vdh_covid19, here::here('data-raw', filename))
  message('Added data from new dates: ', newdates)
}
