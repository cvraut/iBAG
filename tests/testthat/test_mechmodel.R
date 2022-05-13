library(iBAG)

testthat::test_that("test single gene value mech model", {
  demo_data <- iBAG::iBAG_data$new()
  testthat::expect_silent(res <- iBAG::mech.model(demo_data$get.mrna(),demo_data$get.data()))
})

testthat::test_that("test empty.X.constructor", {
  demo_data <- iBAG::iBAG_data$new()
  actual.X <- empty.X.constructor(patients = demo_data$get.patients(),
                                  genes = demo_data$get.genes(),
                                  k = demo_data$get.n_data(),
                                  data.names = demo_data$get.data_names())
  expected.X <- matrix(NA,nrow = 770,ncol = 3*1218)
  row.names(expected.X) <- demo_data$get.patients()
  genes <- demo_data$get.genes()
  colnames(expected.X) <- as.vector(unname(sapply(c("cnv","meth","other"), function(platform){paste(platform,genes,sep = "_")})))
  testthat::expect_identical(actual.X,expected.X,label = "test that empty.X works with demo data")
})
