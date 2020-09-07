#### Tests for filter model scenario functions.

### Create a data frame that defines a Scenario

scenario_test_input <- dplyr::bind_rows(
    ## alpha parameter, two localities and a default, mix of relative and
    ## absolute offsets.
    data.frame(locality='A', time=c(1, 2, 4, 8), isrelative=TRUE,
               value=c(1.2,1.4,1.8,2.0), parm='alpha', stringsAsFactors=FALSE),
    data.frame(locality='B', time=c(1,3,5,7,9), isrelative=FALSE,
               value=c(-1,0,1,2,5), parm='alpha', stringsAsFactors=FALSE),
    data.frame(locality=NA_character_, time=c(2,7), isrelative=TRUE,
               value=c(1.5,1), parm='alpha', stringsAsFactors=FALSE),

    ## beta parameter, one locality, no default, different break times from above
    data.frame(locality='A', time=c(2,4,6,8), isrelative=TRUE,
               value=c(1.5,2.0,1.25,1), parm='beta', stringsAsFactors=FALSE)
)

test_that('scenario construction works', {
  scen <- expect_silent(Scenario(scenario_test_input))
  expect_true(is.Scenario(scen))
  expect_length(scen, 2)
  expect_named(scen, c('alpha', 'beta'), ignore.order=TRUE)
  expect_true(all(sapply(scen, is.data.frame)))


  scen2 <- expect_silent(Scenario(scen))
  expect_equal(scen, scen2)

  expect_error(Scenario(3), regexp='is not TRUE')
  expect_error(Scenario(scenario_test_input[-c(5)]), regexp='is not TRUE')
})


test_that('scenario localization works', {
  scen <- Scenario(scenario_test_input)
  scen_A <- localize_scenario(scen, 'A')
  expect_false(is.Scenario(scen_A))
  expect_s3_class(scen_A, 'LocalScen')
  expect_length(scen_A, 2)
  expect_true(all(sapply(scen, is.data.frame)))

  scen_B <- localize_scenario(scen, 'B')
  expect_false(is.Scenario(scen_B))
  expect_s3_class(scen_B, 'LocalScen')
  expect_length(scen_B, 2)
  expect_s3_class(scen_B[['alpha']], 'data.frame')
  expect_null(scen_B[['beta']])

  scen_C <- localize_scenario(scen, 'C')
  expect_false(is.Scenario(scen_C))
  expect_s3_class(scen_C, 'LocalScen')
  expect_length(scen_C, 2)
  expect_s3_class(scen_C[['alpha']], 'data.frame')
  expect_null(scen_C[['beta']])
})

test_that('scenario parm values are extracted', {
  timevals <- seq(0,10)
  scen <- Scenario(scenario_test_input)
  pb <- c(alpha=1, beta=2)

  parma <- scenario_parm_value(localize_scenario(scen, 'A'), timevals, pb)
  expect_is(parma, 'matrix')
  expect_equal(dim(parma), c(11, 2))
  expect_equal(colnames(parma), c('alpha','beta'))
  expect_equal(parma[ ,'alpha'], c(1, 1.2, 1.4, 1.4, 1.8, 1.8, 1.8, 1.8, 2.0, 2.0, 2.0))
  expect_equal(parma[, 'beta'], 2 * c(1, 1, 1.5, 1.5, 2, 2, 1.25, 1.25,
                                      rep(1,3)))

  parmb <- scenario_parm_value(localize_scenario(scen, 'B'), timevals, pb)
  expect_is(parmb, 'matrix')
  expect_equal(dim(parmb), c(11, 2))
  expect_equal(colnames(parmb), c('alpha','beta'))
  expect_equal(parmb[ ,'alpha'], c(1, 0, 0, 1, 1, 2, 2, 3, 3, 6, 6))
  expect_equal(parmb[ ,'beta'], rep(2, 11))

  parmc <- scenario_parm_value(localize_scenario(scen, 'C'), timevals, pb)
  expect_is(parmc, 'matrix')
  expect_equal(dim(parmc), c(11, 2))
  expect_equal(colnames(parmc), c('alpha','beta'))
  expect_equal(parmc[ ,'alpha'], c(rep(1, 2), rep(1.5, 5), rep(1,4)))
  expect_equal(parmc[ ,'beta'], rep(2, 11))
})

test_that('change times are extracted correctly', {
  scen <- Scenario(scenario_test_input)
  pb <- c(alpha=1, beta=2)
  scenA <- localize_scenario(scen, 'A')
  scenB <- localize_scenario(scen, 'B')
  scenC <- localize_scenario(scen, 'C')

  expect_equal(scenario_change_times(scenA), c(1,2,4,6,8))
  expect_equal(scenario_change_times(scenB), c(1,3,5,7,9))
  expect_equal(scenario_change_times(scenC), c(2, 7))
})

