library(iBAG)

testthat::test_that("Check null constructor for iBAG_data",{
  my.data <- iBAG_data$new()

  #testthat::expect_true(all(dim(my.data$get.meth()) == c(163,176)))
  #testthat::expect_true(all(dim(my.data$get.mrna()) == c(163,48)))
  #testthat::expect_true(all(dim(my.data$get.cnv()) == c(163,482)))
  #testthat::expect_true(all(dim(my.data$get.outcome()) == c(163,1)))
})
