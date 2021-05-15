# Data objects for demo and testing purposes

# TODO Complete Documentation for data associated with package.
# Fill in what stuff means and replace all '???' with defined values
# Also verify that the descriptions are correct & accurate

#' Methylation data for ???
#'
#' A dataset containing the methylation data for 163 individuals for ??? genes across ??? sites.
#' Data obtained from ??? at ???
#'
#' \itemize{
#'   \item ???
#' }
#'
#' @docType data
#' @keywords datasets
#' @name meth
#' @usage data(meth)
#' @format A data frame with ??? rows and ??? variables
#'
#' @export
meth<-read.csv("./data/methylationdata.csv")

#' Gene expression data
#'
#' A dataset containing the mrna expression data for 163 individuals for ??? genes across ??? sites.
#' Data obtained from ??? at ???
#'
#' \itemize{
#'   \item ???
#' }
#'
#' @docType data
#' @keywords datasets
#' @name mrna
#' @usage data(mrna)
#' @format A data frame with 163 rows and ??? variables
#'
#' @export
mrna<-read.csv("./data/mrnadata.csv")

#' Copy Number Variation data
#'
#' A dataset containing the copy number variation data for 163 individuals for ??? genes across ??? sites total
#'
#' \itemize{
#'   \item ???
#' }
#'
#' @docType data
#' @keywords datasets
#' @name cnv
#' @usage data(cnv)
#' @format A data frame with ??? rows and ??? variables
#'
#' @export
cnv<-read.csv("./data/copynumberdata.csv")

#' Survival data
#'
#' Survival data, measured in ???, of 163 individuals suffering from Glioblastoma.
#' Data obtained from ??? at ???.
#'
#' \itemize{
#'   \item X. TCGA patient ID
#'   \item V1. Patient outcomes measured in ???
#' }
#'
#' @docType data
#' @keywords datasets
#' @name dsurv
#' @usage data(dsurv)
#' @format A data frame with 163 rows and 2 variables
#'
#' @export
dsurv<-read.csv("./data/survivaltimes.csv")
