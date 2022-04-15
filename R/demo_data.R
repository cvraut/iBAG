# demo_data.R
# this File generates the demo data
# by default it always loads the demo data

#' Methylation data for ???
#'
#' @docType data
#' @keywords datasets
#' @name demo_meth
#' @usage data(demo_meth)
#' @format A data frame with ??? rows and ??? variables
#'
#' @export
demo_meth <-read.csv(system.file("data","methylationdata.csv",package = "iBAGpkg"), row.names = 1,header = TRUE)

#' mRNA expression data for ???
#'
#' @docType data
#' @keywords datasets
#' @name demo_mrna
#' @usage data(demo_mrna)
#' @format A data frame with ??? rows and ??? variables
#'
#' @export
demo_mrna <-read.csv(system.file("data","mrnadata.csv",package = "iBAGpkg"), row.names = 1,header = TRUE)

#' cnv data for ???
#'
#' @docType data
#' @keywords datasets
#' @name demo_cnv
#' @usage data(demo_cnv)
#' @format A data frame with ??? rows and ??? variables
#'
#' @export
demo_cnv <-read.csv(system.file("data","copynumberdata.csv",package = "iBAGpkg"), row.names = 1,header = TRUE)

#' outcome data for ???
#'
#' @docType data
#' @keywords datasets
#' @name demo_outcome
#' @usage data(demo_outcome)
#' @format A data frame with ??? rows and ??? variables
#'
#' @export
demo_outcome <-read.csv(system.file("data","survivaltimes.csv",package = "iBAGpkg"), row.names = 1,header = TRUE)
