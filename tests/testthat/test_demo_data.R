library(iBAG)

testthat::test_that("Check data sizes & loading",{
  testthat::expect_true(all(dim(iBAG::demo_meth) == c(770,1215)))
  testthat::expect_true(all(dim(iBAG::demo_mrna) == c(770,1215)))
  testthat::expect_true(all(dim(iBAG::demo_cnv) == c(770,1215)))
  testthat::expect_true(all(dim(iBAG::demo_outcome) == c(770,1)))
})

testthat::test_that("Check data columns contain some variance",{
  testthat::expect_true(all(apply(iBAG::demo_cnv,2,sd) != 0))
  testthat::expect_true(length(apply(iBAG::demo_cnv,2,sd) != 0) == 1215)
  testthat::expect_true(all(apply(iBAG::demo_meth,2,sd) != 0))
  testthat::expect_true(all(apply(iBAG::demo_mrna,2,sd) != 0))
  testthat::expect_true(all(apply(iBAG::demo_outcome,2,sd) != 0))
})