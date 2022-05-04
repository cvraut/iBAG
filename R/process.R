# process.R
# file to process & sanitize data

#' data.validate
#' @name data.validate
#' @description Validates mrna, outcome, and upstream data to ensure patients & genes match.
#' @details
#' This function checks:
#'  - patients match across the list of datasets in data and mrna
#'  - patients in outcome are in mrna and data
#'  - genes match across mrna & data
#'
#' A match across patient vector is defined as equivalent order and content.
#' A match across gene vector is defined as each gene in mrna appears at least
#' once in each dataset in data.
#' @usage data.validate(mrna,outcome,data)
#'
#' @param mrna mrna data (follow convention of demo_mrna)
#' @param outcome outcome data (follow demo conventions)
#' @param data list of other dataframes (follow demo conventions)
#' @param DEBUG (FALSE) flag to print debug statements
#' @return (boolean) whether this set of iBAG data is valid
#'
#' @export
data.validate <- function(mrna,
                          outcome,
                          data,
                          DEBUG = FALSE,
                          ...){
  return(FALSE)
}


#' data.validate.patients
#' @name data.validate.patients
#' @description Validates mrna, outcome, and upstream data to ensure patients match.
#' @details
#' This function checks:
#'  - patients match across the list of datasets in data and mrna
#'  - patients in outcome are in mrna and data
#'
#' A match across patient vector is defined as equivalent order and content.
#' @usage data.validate.patients(mrna,outcome,data)
#'
#' For each dataset the patient vector is assumed to be the rownames.
#' @param mrna :mrna data (follow convention of demo_mrna)
#' @param outcome: outcome data (follow demo conventions)
#' @param data :list of other dataframes (follow demo conventions)
#' @param DEBUG FALSE: flag to print debug statements
#' @return boolean: whether this set of iBAG data is valid for the patients only
#'
#'
#' @export
data.validate.patients <- function(mrna,
                                   outcome,
                                   data,
                                   DEBUG = FALSE,
                                   ...){
  pat.list <- row.names(mrna)
  if(!identical(pat.list,row.names(outcome))){
    if(DEBUG){
      print("mrna & outcome ids don't match")
    }
    return(FALSE)
  }
  result <- sapply(data,FUN = function(data){return(identical(pat.list,row.names(data)))})
  if(DEBUG){
    print("output of identical across rows of data:")
    print(result)
  }
  return(all(result))
}

#' data.collapse
#' @name data.collapse
#' @description Collapses data which has >1 column per gene using PCA
#' @usage data.collapse(mrna, data)
#'
#' @param mrna : dataframe for mrna (follows demo)
#' @param data : list of data that needs to consolidated
#' @param PC_VAR_THRESH : The threshold for PCs to accept
#' @param pc_collapse : The additional consolidation function applied to selected PCs
#' @param DEBUG FALSE: flag to print debug statements
#'
#' @export
data.collapse <- function(mrna,
                          data,
                          PC_VAR_THRESH = 0.09,
                          pc_collapse = function(pcs){return(pcs[1])},
                          DEBUG=FALSE,...){
  genes <- colnames(mrna)
  n <- dim(mrna)[1]
  p <- length(genes)

  collapse_gene_data <- function(gene_i,data){
    ind_data <- grep(genes[gene_i],colnames(data))
    if(DEBUG){
      cat(sprintf("probes: %d\n",length(ind_data)))
    }
    if (length(ind_data)==0){
      scores_data <- rep(0,n)
    }  else if (length(ind_data)==1) {
      scores_data <- as.matrix(data[,ind_data])
    } else {
      PCA_data <- princomp(data[,ind_data])
      num_scores_data <- which(cumsum(PCA_data$sdev^2/sum(PCA_data$sdev^2))>=PC_VAR_THRESH)[1]
      if(num_scores_data > 1){
        scores_data  <- as.matrix(apply(PCA_data$scores[,1:num_scores_data],1,pc_collapse))
      } else {
        scores_data  <- as.matrix(PCA_data$scores[,1:num_scores_data],1,pc_collapse)
      }
    }
    if(DEBUG){
      cat(sprintf("post_PCA: %d\n",dim(scores_data)[2]))
    }
    return(scores_data)
  }

  # note parallel later ...
  data.collapsed <- sapply(1:p,FUN=function(i){collapse_gene_data(i,data)})
  row.names(data.collapsed) <- row.names(mrna)
  colnames(data.collapsed) <- colnames(mrna)
  return(data.collapsed)
}
