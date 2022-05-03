# mechmodel.R
# file that runs & fits the mechanistic model

library(mgcv)

#' mechmodel
#'
#' @param meth
#' @param data.list
#' @param DEBUG
#' @param ...
#'
#' @export
mech.model <- function(mrna,data.list,DEBUG=TRUE,...){
  genes <- colnames(mrna)
  n <- dim(mrna)[1]
  p <- length(genes)
  k = length(data.list) + 1

  PC_VAR_THRESH = 0.09

  X <- matrix(NA,nrow=n,ncol=p*k)
  row.names(X) <- rownames(mrna)
  colnames(X) <- paste(rep(sapply(1:k,FUN = function(i){sprintf(ifelse(i<k,"data_%d","other"),i)}),each=p),rep(genes,k),sep="_")
  SS <- matrix(NA,nrow=p,ncol = k+1)
  colnames(SS) <- c("SST",sapply(1:(k-1),FUN = function(i){sprintf("SSD%d",i)}),"SSO")
  row.names(SS) <- genes

  # TODO: check if gene data has been collapsed & collapse it if needed
  get_gene <- function(gene_i,data.list){
    all_data <- sapply(1:(k-1), function(data_i){
      ind_data <- grep(genes[gene_i],colnames(data.list[[data_i]]));
      return(data.list[[data_i]][,ind_data])
    })
    return(all_data)
  }


  process_gene <- function(gene_i){
    if(mean(mrna[,gene_i]) != 0){
      mrna[,gene_i] <- mrna[,gene_i] - mean(mrna[,gene_i])
      if(DEBUG){
        cat("Someone did not zero-mean the mrna >:(\n")
      }
    }

    scores_data <- get_gene(gene_i,data.list)
    scores.form <- paste(sapply(1:(k-1),FUN=function(i){sprintf("s(scores_data[,%d])",i)}),collapse = ' + ')
    mrna.1 <- mrna[,gene_i]
    formula_all <- paste(sprintf("mrna.1 ~ -1 + "),scores.form)
    if(DEBUG){
      cat(sprintf("Gene: %d\nFormula: %s\n",gene_i,formula_all))
    }
    gam.mrna  <- gam(as.formula(formula_all))
    fit_scores <- as.matrix(predict.gam(gam.mrna,type="terms"))
    if(DEBUG){
      cat("gam_fit done!")
    }
    sapply(1:(k-1),FUN = function(i){
      X[,paste(sprintf("data_%d",i),genes[gene_i],sep="_")] <<- fit_scores[,i];
      SS[gene_i,i+1] <<- sum( ( (fit_scores[,i]) - mean(mrna[,gene_i]) )^2  )
    })

    O <- gam.mrna$residuals
    X[,paste("other",genes[gene_i],sep="_")] <<- O
    # Pseudo Sums of Squares (to use to find percentages of explained variance)
    SS[gene_i,1] <<- sum( (mrna[,gene_i] - mean(mrna[,gene_i]))^2 )
    if(k>2){
      SS[gene_i,(k+1)] <<- SS[gene_i,1] - rowSums(SS[gene_i,2:k])
    } else {
      SS[gene_i,3] <<- SS[gene_i,1] - SS[gene_i,2]
    }
    if(DEBUG){
      print(dim(SS))
    }
  }
  # note parallel later ...
  sapply(1:p,FUN=function(i){process_gene(i)})
  return(list(X=X,
              SS=SS))
}
