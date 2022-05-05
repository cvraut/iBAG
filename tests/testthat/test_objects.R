library(iBAG)

testthat::test_that("Check null constructor for iBAG_data",{
  my.data <- iBAG_data$new()

  testthat::expect_true(all(dim(my.data$get.data("meth")) == c(770,1218)))
  testthat::expect_true(all(dim(my.data$get.mrna()) == c(770,1218)))
  testthat::expect_true(all(dim(my.data$get.data("cnv")) == c(770,1218)))
  testthat::expect_true(all(dim(my.data$get.outcome()) == c(770,1)))
  testthat::expect_true(my.data$get.n_genes() == 1218)
  testthat::expect_true(my.data$get.n_patients() == 770)
  testthat::expect_true(my.data$get.n_genes() == 1218)
  testthat::expect_true(my.data$get.n_data() == 2)
})
