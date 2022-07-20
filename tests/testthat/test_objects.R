library(iBAG)

testthat::test_that("Check null constructor for iBAG_data (getters & setters)", { # nolint
  my.data <- iBAG_data$new()

  testthat::expect_true(all(dim(my.data$get.data("meth")) == c(770,  1215)))
  testthat::expect_true(all(dim(my.data$get.mrna()) == c(770, 1215)))
  testthat::expect_true(all(dim(my.data$get.data("cnv")) == c(770, 1215)))
  testthat::expect_true(all(dim(my.data$get.outcome()) == c(770, 1)))
  testthat::expect_true(my.data$get.n_genes() == 1215)
  testthat::expect_true(my.data$get.n_patients() == 770)
  testthat::expect_true(my.data$get.n_genes() == 1215)
  testthat::expect_true(my.data$get.n_data() == 2)
})

testthat::test_that("Check construction of empty iBAG_results", {
  testthat::expect_silent(my.results <- iBAG_results$new())

  testthat::expect_identical(my.results$get.X(), matrix(0, 1, 1))
  testthat::expect_identical(my.results$get.Y(), c(0))
  testthat::expect_identical(my.results$get.SS(), list(c(0), c(0), c(0)))
  testthat::expect_identical(my.results$get.beta_mean(), c(0))
  testthat::expect_identical(my.results$get.beta_incl_prob(), c(0))
  testthat::expect_equal(my.results$get.n_patients(), 1)
  testthat::expect_equal(my.results$get.n_data(), 1)
  testthat::expect_equal(my.results$get.n_genes(), 1)
})

testthat::test_that("Check getters of iBAG_results", {
  my.results <- iBAG_results$new(X = matrix("X", 1, 1),
                                 Y = 2,
                                 SS = list(SST=4.4, 345, SSO="zero"),
                                 beta_mean="b_mean",
                                 beta_incl_prob=TRUE,
                                 validate = FALSE)

  testthat::expect_equal(my.results$get.X(), matrix("X", 1, 1))
  testthat::expect_equal(my.results$get.Y(), 2)
  testthat::expect_identical(my.results$get.SS(), list(SST=4.4, 345, SSO="zero"))
  testthat::expect_equal(my.results$get.SS(1), 4.4)
  testthat::expect_equal(my.results$get.SS("SSO"), "zero")
  testthat::expect_equal(my.results$get.beta_mean(), "b_mean")
  testthat::expect_equal(my.results$get.beta_incl_prob(), TRUE)
  testthat::expect_equal(my.results$get.n_patients(), 1)
  testthat::expect_equal(my.results$get.n_data(), 1)
  testthat::expect_equal(my.results$get.n_genes(), 1)
})

testthat::test_that("Check validate of iBAG_results", {
  testthat::expect_warning(iBAG_results$new(X = NULL), "X, Y, or SS was missing. Setting validate to FALSE")
  testthat::expect_error(iBAG_results$new(X=matrix(1, nrow=2, ncol=1)), "n_patients is not consistent across X and Y")
  testthat::expect_error(iBAG_results$new(SS=list(c(1, 2), c(0))), "SS does not contain at least 3 elements")
  testthat::expect_error(iBAG_results$new(SS=list(c(1, 2), c(0), c(0), c(0))),
    "X does not contain the right amount of gene columns according to SS")
  testthat::expect_error(iBAG_results$new(SS=list(c(1, 2), c(0), c(0))),
    "SS does not contain the correct amount of data for each gene")
  testthat::expect_error(iBAG_results$new(beta_mean = 1:2), "dimensions of beta_mean not consistent with X")
  testthat::expect_error(iBAG_results$new(beta_incl_prob = 1:2), "dimensions of beta_incl_prob not consistent with X")
})

testthat::test_that("Test get.data.names from iBAG_data",{
  testthat::expect_equal(iBAG_data$new()$get.data_names(),c("cnv","meth"),
                         label = "make sure demo data works")
})
