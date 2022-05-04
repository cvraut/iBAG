library(iBAG)

testthat::test_that("Check data sizes & loading",{
  testthat::expect_true(all(dim(demo_meth) == c(770,1218)))
  testthat::expect_true(all(dim(demo_mrna) == c(770,1218)))
  testthat::expect_true(all(dim(demo_cnv) == c(770,1218)))
  testthat::expect_true(all(dim(demo_outcome) == c(770,1)))
})
