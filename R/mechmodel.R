# mechmodel.R
# file that runs & fits the mechanistic model

library(mgcv)

#' mechmodel
#'
#' @param meth
#' @param mrna
#' @param cnv
#' @param DEBUG
#' @param ...
#'
#' @export
mechmodel <- function(meth,mrna,cnv,DEBUG=TRUE,...){
  genes <- colnames(mrna)
  n <- dim(mrna)[1]
  p <- length(genes)
  k = 3

  PC_VAR_THRESH = 0.09

  X <- matrix(NA,nrow=n,ncol=p*k)
  row.names(X) <- rownames(mrna)
  colnames(X) <- paste(rep(c("Meth","CN","Other"),each=p),rep(genes,3),sep="_")
  SS <- matrix(NA,nrow=p,ncol = 4)
  colnames(SS) <- c("SST","SSM","SSCN","SSO")
  row.names(SS) <- genes
  collapse_gene_data <- function(gene_i,data){
    ind_data <- grep(genes[gene_i],colnames(data))
    if(DEBUG){
      cat(sprintf("probes: %d\n",length(ind_data)))
    }
    if (length(ind_data)==0){
      scores_data <- rep(0,n)
    }  else if (length(ind_data)==1) {
      scores_data <- as.matrix(data[,ind_data])
    } else {  ## If only 1 data value, keep raw data (no PCA).
        PCA_data <- princomp(data[,ind_data])
        num_scores_data <- which(cumsum(PCA_data$sdev^2/sum(PCA_data$sdev^2))>=PC_VAR_THRESH)[1]
        scores_data  <- as.matrix(PCA_data$scores[,1:num_scores_data])
    }
    if(DEBUG){
      cat(sprintf("post_PCA: %d\n",dim(scores_data)[2]))
    }
    return(scores_data)
  }
  process_gene <- function(gene_i){
    if(mean(mrna[,gene_i]) != 0){
      mrna[,gene_i] <- mrna[,gene_i] - mean(mrna[,gene_i])
      if(DEBUG){
        cat("Someone did not zero-mean the mrna >:(\n")
      }
    }
    scores_meth <- collapse_gene_data(gene_i,meth)
    num_scores_meth <- dim(scores_meth)[2]
    scores_cnv <- collapse_gene_data(gene_i,meth)
    num_scores_cnv <- dim(scores_cnv)[2]

    if (num_scores_meth == 1){
      formula_meth <- "s(scores_meth[,1])"
    }  else{
      formula_meth <- paste("s(scores_meth[,",paste(1:num_scores_meth, collapse="]) + s(scores_meth[,"),"])",sep='')
    }
    if (length(scores_cnv) == n){
      formula_cnv <- "s(scores_cnv[,1])"
    } else{
      formula_cnv <- paste("s(scores_cnv[,",paste(1:num_scores_cnv, collapse="]) + s(scores_cnv[,"),"])",sep='')
    }
    mrna.1 <- mrna[,gene_i]
    formula_all <- paste("mrna.1 ~ -1 + ",formula_meth," + ",formula_cnv)
    if(DEBUG){
      cat(sprintf("Gene: %d\nFormula: %s\n",gene_i,formula_all))
    }
    gam.mrna  <- gam(as.formula(formula_all))
    fit_meth <- as.matrix(predict.gam(gam.mrna,type="terms")[,1:num_scores_meth])
    fit_cnv <- as.matrix(predict.gam(gam.mrna,type="terms")[,(num_scores_meth+1):(num_scores_meth+num_scores_cnv)])

    M <- apply(fit_meth,1,sum)
    CN <- apply(fit_cnv,1,sum)
    O <- gam.mrna$residuals
    X[,paste("Meth",genes[gene_i],sep="_")]   <-  M
    X[,paste("CN",genes[gene_i],sep="_")] <- CN
    X[,paste("Other",genes[gene_i],sep="_")] <- O
    # Pseudo Sums of Squares (to use to find percentages of explained variance)
    SS[gene_i,1] <- sum( (mrna[,gene_i] - mean(mrna[,gene_i]))^2 )
    SS[gene_i,2] <- sum( ( (M) - mean(mrna[,gene_i]) )^2  )
    SS[gene_i,3] <- sum( ( (CN) - mean(mrna[,gene_i]) )^2  )
    SS[gene_i,4] <- SS[gene_i,1] - SS[gene_i,2] - SS[gene_i,3]
    if(DEBUG){
      print(SS[gene_i,1])
      print(SS[gene_i,2])
      print(SS[gene_i,3])
      print(SS[gene_i,4])
      if(SS[gene_i,4] < 0){
        print("we got a problem here")
      }
    }
  }
  # note parallel later ...
  sapply(1:p,FUN=process_gene)
  return(list(X=X,
              SS=SS))
}
