# mechmodel.R
# file that runs & fits the mechanistic model

library(mgcv)

#' empty.X.constructor
#'
#' @description This is a helper function to create & name the rows & columns of X appropriately
#' @param k this is the number of datasets (not including the other/residual, gets added automatically)
empty.X.constructor <- function(patients,
                                genes,
                                k,
                                data.names,
                                other.name = "other",
                                sep = "_",
                                DEBUG = FALSE){
  n = length(patients)
  p = length(genes)
  X <- matrix(NA, nrow=n, ncol=p*(k+1))
  row.names(X) <- patients

  if(other.name %in% data.names){
    stop("other.name is identical to one of the dataset names.")
  }

  X.namer <- function(i){
    return(sprintf(ifelse(i<=k,data.names[i],other.name),i))
  }
  colnames(X) <- paste(rep(sapply(1:(k+1),FUN = X.namer),each=p),rep(genes,k),sep="_")
  return(X)
}

#' mech.model
#'
#' @description This fits the mechanistic model for iBAG
#'
#' @param mrna the mrna data. Look at iBAG::demo_mrna for example format
#' @param data.list a list of the upstream data
#' @param sep (_): the character value to use as a seperator
#' @param other.name ("other"): the string to use to label residuals from mechmodel
#' @param default.data.name ("data"): the string to use to label unnamed datasets from data.list
#' @param DEBUG (FALSE): debug flag
#' @param ...
#' @export
mech.model <- function(mrna,
                       data.list,
                       sep = "_",
                       other.name = "other",
                       default.data.name = "data",
                       DEBUG=FALSE, ...){
  genes <- colnames(mrna)
  n <- nrow(mrna)
  p <- length(genes)
  k <- length(data.list) + 1
  data_names <- get.data.names(data.size = k-1,
                               data.names = names(data.list),
                               sep = sep,
                               default.data.name = default.data.name,
                               DEBUG = DEBUG)
  X <- empty.X.constructor(patients = row.names(mrna),
                           genes = genes,
                           k = (k-1),
                           data.names = data_names,
                           other.name = other.name,
                           sep = sep,
                           DEBUG = DEBUG)
  SS <- matrix(NA,nrow=p,ncol = k+1)
  colnames(SS) <- c("SST",sapply(1:(k-1),FUN = function(i){sprintf("SS%s%s",sep,data_names[i])}),"SSO")
  row.names(SS) <- genes

  # TODO: check if gene data has been collapsed & collapse it if needed
  # TODO: support multiple prob information
  get_gene <- function(gene_i,data.list){
    all_data <- sapply(1:(k-1), function(data_i){
      gene_regex <- paste("^",sprintf("%s(%s.+)?",genes[gene_i],sep),"$",sep = "")
      ind_data <- grep(gene_regex,colnames(data.list[[data_i]]));
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
      X[,gene_i+(i-1)*p] <<- fit_scores[,i];
      SS[gene_i,i+1] <<- sum( ( (fit_scores[,i]) - mean(mrna[,gene_i]) )^2  )
    })

    O <- gam.mrna$residuals
    X[,gene_i+(k-1)*p] <<- O
    # Pseudo Sums of Squares (to use to find percentages of explained variance)
    SS[gene_i,1] <<- sum( (mrna[,gene_i] - mean(mrna[,gene_i]))^2 )
    if(k>2){
      SS[gene_i,(k+1)] <<- SS[gene_i,1] - sum(SS[gene_i,2:k])
    } else {
      SS[gene_i,3] <<- SS[gene_i,1] - SS[gene_i,2]
    }
    if(DEBUG){
      print(dim(SS))
    }
  }
  # TODO: parallel later ...
  sapply(1:p,FUN=function(i){process_gene(i)})
  return(list(X=X,
              SS=SS))
}
