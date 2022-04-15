library(iBAGpkg)

testthat::test_that("Check data sizes & loading",{
  testthat::expect_true(all(dim(demo_meth) == c(163,176)))
  testthat::expect_true(all(dim(demo_mrna) == c(163,48)))
  testthat::expect_true(all(dim(demo_cnv) == c(163,482)))
  testthat::expect_true(all(dim(demo_outcome) == c(163,1)))
})
