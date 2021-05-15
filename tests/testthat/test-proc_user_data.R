context("proc_user_data")

testthat::test_that("Null Checks: Errors",{
  testthat::expect_error(proc_user_data(NULL,NULL,NULL,NULL),
                         "Missing gene expression data")


})
