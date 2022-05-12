library(iBAG)

testthat::test_that("test single gene value mech model", {
  demo_data <- iBAG::iBAG_data$new()
  testthat::expect_silent(res <- iBAG::mech.model.single_gene_val(demo_data$get.mrna(),demo_data$get.data()))
})

testthat::test_that("test empty.X.constructor", {
  demo_data <- iBAG::iBAG_data$new()
  patients = demo_data$get.patients(),
  genes = demo_data$get.genes(),
  k = demo_data$get.n_data(),
  data.names = demo_data$get.data_names()
})
