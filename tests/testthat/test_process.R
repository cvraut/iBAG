library(iBAG)

testthat::test_that("Check data.validate.patients",{
  test.patients.0 = sapply(1:20,FUN=function(i){sprintf("pat_%d",i)})
  test.patients.1 = sapply(1:20,FUN=function(i){sprintf("pat_%d",i)})
  test.patients.2 = sapply(0:19,FUN=function(i){sprintf("pat_%d",i)})
  test.patients.3 = sapply(1:21,FUN=function(i){sprintf("pat_%d",i)})
  test.patients.4 = sapply(1:19,FUN=function(i){sprintf("pat_%d",i)})
  test.patients.5 = sapply(20:1,FUN=function(i){sprintf("pat_%d",i)})

  #datas
  data1 <- matrix(NA,nrow = 20,ncol = 1)
  data2 <- matrix(NA,nrow = 20,ncol = 1)

  testthat::expect_true(data.validate.patients(mrna=data1,outcome=data2,data=list(data1,data2)),
                        label = "normal use (rownames are integers), everything matches")

  row.names(data1) <- test.patients.0
  row.names(data2) <- test.patients.1
  data3 <- matrix(NA,nrow = 20,ncol = 1)
  row.names(data3) <- test.patients.2
  data4 <- matrix(NA,nrow = 21,ncol = 1)
  row.names(data4) <- test.patients.3
  data5 <- matrix(NA,nrow = 19,ncol = 1)
  row.names(data5) <- test.patients.4
  data6 <- matrix(NA,nrow = 20,ncol = 1)
  row.names(data6) <- test.patients.5

  testthat::expect_true(data.validate.patients(mrna=data1,outcome=data2,data=list(data1,data2)),
                        label = "normal use (rownames are strings), everything matches")

  #bad cases
  testthat::expect_false(data.validate.patients(mrna=data1,outcome=data2,data=list(data1,data3)),
                         label = "(rownames are strings) one dataset is off but same length")
  testthat::expect_false(data.validate.patients(mrna=data1,outcome=data2,data=list(data4,data2)),
                         label = "(rownames are strings) one dataset is off but longer")
  testthat::expect_false(data.validate.patients(mrna=data1,outcome=data5,data=list(data4,data2)),
                         label = "(rownames are strings) outcome is off but shorter")
  testthat::expect_false(data.validate.patients(mrna=data6,outcome=data2,data=list(data1,data2)),
                         label = "(rownames are strings) mrna is off but reversed")
})
