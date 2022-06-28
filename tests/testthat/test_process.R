library(iBAG)

testthat::test_that("Check data.validate.patients",{
  test.patients.0 = sapply(1:20,FUN=function(i){sprintf("pat_%d",i)})
  test.patients.1 = sapply(1:20,FUN=function(i){sprintf("pat_%d",i)})
  test.patients.2 = sapply(0:19,FUN=function(i){sprintf("pat_%d",i)})
  test.patients.3 = sapply(1:21,FUN=function(i){sprintf("pat_%d",i)})
  test.patients.4 = sapply(1:19,FUN=function(i){sprintf("pat_%d",i)})
  test.patients.5 = sapply(20:1,FUN=function(i){sprintf("pat_%d",i)})
  test.patients.6 = sapply(c(1:10,1:10),FUN=function(i){sprintf("pat_%d",i)})

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
  data7 <- matrix(NA,nrow = 20,ncol = 1)
  row.names(data7) <- test.patients.6

  testthat::expect_true(data.validate.patients(mrna=data1,outcome=data2,data=list(data1)),
                        label = "normal use (rownames are strings), everything matches - 1 dataset")

  testthat::expect_true(data.validate.patients(mrna=data1,outcome=data2,data=list(data1,data2)),
                        label = "normal use (rownames are strings), everything matches - 2 datasets")

  #bad cases
  testthat::expect_false(data.validate.patients(mrna=data1,outcome=data2,data=list(data3)),
                         label = "(rownames are strings) one dataset is off but same length - 1 dataset")
  testthat::expect_false(data.validate.patients(mrna=data1,outcome=data2,data=list(data1,data3)),
                         label = "(rownames are strings) one dataset is off but same length - 2 datasets")
  testthat::expect_false(data.validate.patients(mrna=data1,outcome=data2,data=list(data4,data2)),
                         label = "(rownames are strings) one dataset is off but longer")
  testthat::expect_false(data.validate.patients(mrna=data1,outcome=data5,data=list(data4,data2)),
                         label = "(rownames are strings) outcome is off but shorter")
  testthat::expect_false(data.validate.patients(mrna=data6,outcome=data2,data=list(data1,data2)),
                         label = "(rownames are strings) mrna is off but reversed")
  testthat::expect_false(data.validate.patients(mrna=data7,outcome=data7,data=list(data7,data7)),
                         label = "(rownames are strings) we have duplicate patients in all datasets")
  # demo data should pass
  testthat::expect_true(data.validate.patients(mrna=demo_mrna,outcome=demo_outcome,data=list(meth=demo_meth,cnv=demo_cnv)),
                        label = "demo dataset should pass :|")
})


testthat::test_that("Check data.validate.genes",{
  set.seed(1234)
  genes.list <- c(1:5)
  probes.list <- c(1:5)
  sep = "_"
  data.genes.list <- unlist(sapply(genes.list,FUN=function(gene){return(sapply(sample(probes.list,sample(seq_len(length(probes.list)),1),replace = F),FUN = function(probe){paste(gene,sep,probe,sep = "")}))}))

  data1 <- matrix(NA,nrow=2,ncol=5)
  colnames(data1) <- genes.list
  data2 <- matrix(NA,nrow=2,ncol=5)
  colnames(data2) <- genes.list
  data3 <- matrix(NA,nrow=2,ncol=5)
  colnames(data3) <- c(1,2,3,4,1)
  data4 <- matrix(NA,nrow=2,ncol=6)
  colnames(data4) <- c(1,2,3,4,5,1)
  testthat::expect_true(data.validate.genes(mrna=data1,data=list(meth=data1,cnv=data2),sep=NULL),
                        label = "sep = NULL perfect matches - 2 datasets")
  testthat::expect_true(data.validate.genes(mrna=data1,data=list(data1),sep=NULL),
                        label = "sep = NULL perfect matches - 1 datasets")
  testthat::expect_false(data.validate.genes(mrna=data3,data=list(data1,data2),sep=NULL),
                        label = "repeat gene in mrna - 2 datasets")
  testthat::expect_false(data.validate.genes(mrna=data1,data=list(meth=data4),sep=NULL),
                         label = "repeat gene in data - 1 datasets")
  #TODO: finish adding other non-null testcases later
  testthat::expect_true(data.validate.genes(mrna=demo_mrna,data=list(meth=demo_meth,cnv=demo_cnv),sep="_"),
                        label = "demo dataset should pass :|")
})

testthat::test_that("Verify the demo data is valid",{
  testthat::expect_true(data.validate(mrna=demo_mrna,outcome = demo_outcome,data=list(meth=demo_meth,cnv=demo_cnv),sep=NULL),
                        label = "demo dataset should pass :|")
})

testthat::test_that("test get.data.names", {
  # data.size,data.names,sep,default.data.name
  testthat::expect_identical(get.data.names(4,NULL), c("data_1","data_2","data_3","data_4"),
                             "testing if data.list has no names")
  testthat::expect_identical(get.data.names(2,names(iBAG::iBAG_data$new()$get.data())), c("cnv","meth"),
                             "testing the demo data")
  testthat::expect_identical(get.data.names(7,c("my.data","my.data","","","meth","meth","cnv"),sep=""),
                             c("my.data1","my.data2","data1","data2","meth1","meth2","cnv"),
                             "testing a far more complicated example")
  testthat::expect_error(get.data.names(3,c("data","data1","data"),sep=""),
                         "collision on data renaming 'data1' Check sep.")
})

testthat::test_that("test data.validate.columns", {
  # good data
  data1 <- matrix(seq(1,10),nrow = 2, ncol = 5)
  data2 <- matrix(seq(1,8),nrow = 8,ncol = 1)
  # bad data
  data3 <- matrix(c(1,2,3,3),nrow = 2, ncol = 2)
  data4 <- matrix(c(1,1,3,4),nrow = 2, ncol = 2)
  data5 <- matrix(c(1,2),nrow=1,ncol=2)

  testthat::expect_true(iBAG::data.validate.columns(data1,list(cnv=data2)),
    label = "no problems")
  testthat::expect_false(iBAG::data.validate.columns(data1,list(cnv=data2,meth=data3)),
    label = "problem is in data[[2]]")
  testthat::expect_false(iBAG::data.validate.columns(data1,list(cnv=data4,meth=data2)),
    label = "problem is in data[[1]]")
  testthat::expect_false(iBAG::data.validate.columns(data5,list(cnv=data2,meth=data2)),
    label = "problem in mrna")
})