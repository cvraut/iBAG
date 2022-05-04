# demo_data.R
# this File generates the demo data
# by default it always loads the demo data

#' demo_meth
#' @details
#' Methylation data for ???
#'
#' @docType data
#' @keywords datasets
#' @name demo_meth
#' @usage data(demo_meth)
#' @format A data frame with 770 rows and 1218 variables
#'
#' @export
demo_meth <-read.csv(system.file("data","meth.csv",package = "iBAG"),row.names = 1,header = TRUE)

#' demo_mrna
#' @details
#' mRNA expression data for TCGA-BRCA, note we filter out to only 1218 genes to keep the dataset "small".
#' The 1st column refers to the patient ids
#'
#' @docType data
#' @keywords datasets
#' @name demo_mrna
#' @usage data(demo_mrna)
#' @format A data frame with 770 rows and 1218 variables
#'
#' @export
demo_mrna <-read.csv(system.file("data","mrna.csv",package = "iBAG"),row.names = 1,header = TRUE)

#' demo_cnv
#' @details
#' cnv data for TCGA-BRCA, note we filter out to only 1218 genes to keep the dataset "small"
#' We already condense the cnv to single number summaries per gene
#' The 1st column refers to the patient ids
#'
#' @docType data
#' @keywords datasets
#' @name demo_cnv
#' @usage data(demo_cnv)
#' @format A data frame with 770 rows and 1218 variables
#'
#' @export
demo_cnv <-read.csv(system.file("data","cnv.csv",package = "iBAG"),row.names = 1,header = TRUE)

#' demo_outcome
#' @details
#' EREG-mRNA stemness index outcome data for TCGA-BRCA
#'
#' @docType data
#' @keywords datasets
#' @name demo_outcome
#' @usage data(demo_outcome)
#' @format A data frame with 770 rows and 1 variables
#'
#' @export
demo_outcome <-read.csv(system.file("data","outcome.csv",package = "iBAG"),row.names = 1,header = TRUE)
