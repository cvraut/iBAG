library(iBAG)

testthat::test_that("test spike & slab priors", {
  # default
  expected <- list(r_tau = 0.001,
                   d_tau = 0.001,
                   d_pi = 1,
                   r_pi = 1,
                   r_gamma = 0.001,
                   d_gamma = 0.001)

  testthat::expect_equal(iBAG::construct.priors(),expected, 
    label = "default should be the uniformative priors from graper paper")
  
  # setting priors
  expected <- list(r_tau = 1,
                   d_tau = 2,
                   d_pi = 3,
                   r_pi = 4,
                   r_gamma = 5,
                   d_gamma = 6)

  testthat::expect_equal(iBAG::construct.priors(1,2,3,4,5,6),expected, 
    label = "testing setting priors")

  # test bad prior

  testthat::expect_error(iBAG::construct.priors(1,2,3,4,5,0),
                         regexp = "All prior params must be > 0. supplied params: .+(, .+)*")
  testthat::expect_error(iBAG::construct.priors(1,2,-100,4,5,1),
                         regexp = "All prior params must be > 0. supplied params: .+(, .+)*")
})

