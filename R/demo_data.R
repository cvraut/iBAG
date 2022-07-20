# demo_data.R
# this File generates the demo data
# by default it always loads the demo data
# TODO: run performance checks to see if we should always load the test data.

#' demo_meth
#' @details
#' Methylation data for TCGA-BRCA, note we filter out to only 1218 genes to keep the dataset "small".
#' rownames are the patient ids.
#' colnames are the gene names.
#' We also collapse the multiple probe information from TCGA using dataset.collapse.pc.singular so we only get 1 column per gene.
#'
#' @docType data
#' @keywords datasets
#' @name demo_meth
#' @usage data(demo_meth)
#' @format A data frame with 770 rows and 1218 variables
#' @source https://portal.gdc.cancer.gov/
demo_meth <-read.csv(system.file("data","meth.csv",package = "iBAG"),row.names = 1,header = TRUE)

#' demo_mrna
#' @details
#' mRNA expression data for TCGA-BRCA, note we filter out to only 1218 genes to keep the dataset "small".
#' rownames are the patient ids.
#' colnames are the gene names.
#'
#' @docType data
#' @keywords datasets
#' @name demo_mrna
#' @usage data(demo_mrna)
#' @format A data frame with 770 rows and 1218 variables
#' @source https://portal.gdc.cancer.gov/
demo_mrna <-read.csv(system.file("data","mrna.csv",package = "iBAG"),row.names = 1,header = TRUE)

#' demo_cnv
#' @details
#' cnv data for TCGA-BRCA, note we filter out to only 1218 genes to keep the dataset "small"
#' We already condense the cnv to single number summaries per gene from TCGA.
#' rownames are the patient ids.
#' colnames are the gene names.
#'
#' @docType data
#' @keywords datasets
#' @name demo_cnv
#' @usage data(demo_cnv)
#' @format A data frame with 770 rows and 1218 variables
#' @source https://portal.gdc.cancer.gov/
demo_cnv <-read.csv(system.file("data","cnv.csv",package = "iBAG"),row.names = 1,header = TRUE)

#' demo_outcome
#' @details
#' EREG-mRNA stemness index outcome data for TCGA-BRCA. Consult https://doi.org/10.1016/j.cell.2018.03.034 for more details.
#'
#' @docType data
#' @keywords datasets
#' @name demo_outcome
#' @usage data(demo_outcome)
#' @format A data frame with 770 rows and 1 variables
#' @source https://portal.gdc.cancer.gov/
demo_outcome <-read.csv(system.file("data","outcome.csv",package = "iBAG"),row.names = 1,header = TRUE)
