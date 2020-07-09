test_that('mask_indicator works for all types', {
  td <- as.Date(c('2020-04-30', '2020-05-28', '2020-05-29', '2020-06-30'))
  tdiff <- td - as.Date('2020-01-01')
  tn <- as.numeric(tdiff)
  ti <- as.integer(tn)
  
  ans <- c(0,0,1,1)
  expect_equal(mask_indicator(td), ans)
  expect_equal(mask_indicator(tdiff), ans)
  expect_equal(mask_indicator(tn), ans)
  expect_equal(mask_indicator(ti), ans)
})
