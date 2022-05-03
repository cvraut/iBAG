# process.R
# file to process & sanitize data

#' data.collapse
#'
#' @param mrna
#' @param data
#' @param DEBUG
#' @param ...
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
