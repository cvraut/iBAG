# process.R
# file to process & sanitize data

#' data.validate
#' @name data.validate
#' @description Validates mrna, outcome, and upstream data to ensure patients & genes match.
#' @details
#' Runs the data.validate.patients & data.validate.genes subroutines. Consult them directly for more information
#'
#' @usage data.validate(mrna,outcome,data)
#'
#' @param mrna mrna data (follow convention of demo_mrna)
#' @param outcome outcome data (follow demo conventions)
#' @param data list of other dataframes (follow demo conventions)
#' @param sep ("_"): the character seperator to use to differentiate mutliple columns of data per gene
#' @param DEBUG (FALSE) flag to print debug statements
#' @return (boolean) whether this set of iBAG data is valid
#'
#' @export
data.validate <- function(mrna,
                          outcome,
                          data,
                          sep="",
                          DEBUG = FALSE,
                          ...){
  return(iBAG::data.validate.patients(mrna = mrna,outcome = outcome,data = data,DEBUG = DEBUG) &&
         iBAG::data.validate.genes(mrna = mrna,data = data,sep=sep,DEBUG=DEBUG) &&
         iBAG::data.validate.columns(mrna = mrna,data = data,DEBUG=DEBUG) &&
         length(data) > 0)
}


#' data.validate.patients
#' @name data.validate.patients
#' @description Validates mrna, outcome, and upstream data to ensure patients match.
#' @details
#' This function checks:
#'   - patients are unique in mrna
#'   - patients match across the list of datasets in data and mrna
#'   - patients in outcome are in mrna and data
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
  if(length(pat.list) != length(unique(pat.list))){
    if(DEBUG){
      print("mrna has duplicate rows")
    }
    return(FALSE)
  } else if(!identical(pat.list,row.names(outcome))){
    if(DEBUG){
      print("mrna & outcome ids don't match")
    }
    return(FALSE)
  }
  result <- sapply(data,FUN = function(dataset){return(identical(pat.list,row.names(dataset)))})
  if(DEBUG){
    print("output of identical across rows of data:")
    print(result)
  }
  return(all(result))
}

#' data.validate.genes
#' @name data.validate.genes
#' @description Validates mrna and upstream data to ensure genes match.
#' @details
#' This function checks:
#'   - all genes are unique in mrna
#'   - all columns are unique across all datasets
#'   - given an upstream data platform:
#'     - for each gene in mrna there is at least 1 gene present
#'     - upstream platforms do not contain genes which are not present in mrna
#'
#' A match across gene makes use of the {sep} argument. First a regex is constructed
#' from the mrna such as r"{GENE}({SEP}.+)?" then all gene data is found by matching this pattern.
#'
#' Note: gene names are automatically converted to strings
#'
#' Note: SEP should not be able to be found in any of the genes from mrna (auto-fail if sep="")
#'
#' Note: if SEP == NULL then this method will perform perfect matches similar to data.validate.patient
#'
#' Note: the order of the genes don't have to match, but it is important that the contents do.
#'
#' @usage data.validate.genes(mrna,data)
#'
#' For each dataset the patient vector is assumed to be the rownames.
#' @param mrna :mrna data (follow convention of demo_mrna)
#' @param data :list of other dataframes (follow demo conventions)
#' @param sep ("_"): the character seperator to use to differentiate mutliple columns of data per gene
#' @param DEBUG FALSE: flag to print debug statements
#' @return boolean: whether this set of iBAG data is valid for the genes only
#'
#'
#' @export
data.validate.genes <- function(mrna,
                                data,
                                sep="_",
                                DEBUG = FALSE,
                                ...){
  genes.list <- paste(colnames(mrna))
  if(length(genes.list) != length(unique(genes.list)) | any(sapply(data,FUN = function(src){return(length(colnames(src))!=length(unique(colnames(src))))}))){
    if(DEBUG){
      print("genes in mrna or columns of data were not unique")
    }
    return(FALSE)
  }
  if(is.null(sep)){
    genes.list <- sort(genes.list)
    result <- sapply(data,FUN = function(data){return(identical(genes.list,sort(paste(colnames(data)))))})
    if(DEBUG){
      print("output of identical across rows of data:")
      print(result)
    }
    return(all(result))
  } else if(sep == ""){
    if(DEBUG){
      print('sep == "" case')
    }
    return(FALSE)
  } else{
    # if sep is found in genes.list
    if(any(grep(sep,genes.list,fixed=TRUE))){
      return(FALSE)
    }
    get_gene_cnts <- function(gene,sep,probes){
      patt <- sprintf("^%s(%s.+)?$",gene,sep)
      res <- stringr::str_count(probes,pattern = patt)
      return(sum(res))
    }
    check_dataset <- function(dataset){
      probes <- paste(colnames(dataset))
      return(sum(sapply(genes.list,FUN = function(gene){return(get_gene_cnts(gene,sep,probes))})) == length(probes))
    }
    return(all(sapply(data,FUN = function(dataset){return(check_dataset(dataset))})))
  }
}


#' data.validate.columns
#' @name data.validate.columns
#' @description Validates the columns of mrna and upstream data to ensure all columns have sd > 0.
#' @details
#' This function checks:
#'   - all columns in each dataset contain some variation
#'   - returns false if it encounters a column in any dataset without variation
#'     - otherwise returns true
#' 
#' @usage data.validate.columns(mrna,data)
#'
#' @param mrna :mrna data (follow convention of demo_mrna)
#' @param data :list of other dataframes (follow demo conventions)
#' @param DEBUG FALSE: flag to print debug statements
#' @return boolean: whether this set of iBAG data is valid for the columns only
#'
#'
#' @export
data.validate.columns <- function(mrna,
                                  data,
                                  DEBUG = FALSE,
                                  ...){
  check_dataset <- function(data){
    col_sds <- apply(data,2,sd)
    return(all(is.na(col_sds) == FALSE) && all(col_sds != 0))
  }
  if(!check_dataset(mrna)){
    if(DEBUG){
      print("mrna has a column with 0 sd")
    }
    return(FALSE)
  }
  if(any(sapply(data,FUN=check_dataset) == FALSE)){
    if(DEBUG){
      print("upstream data has a column with 0 sd")
    }
    return(FALSE)
  }
  return(TRUE)
}

#' dataset.collapse.pc.singular
#' @name dataset.collapse.pc.singular
#' @description Collapses data which has >1 column per gene using PCA
#' @usage data.collapse(mrna, data)
#'
#' @param mrna : dataframe for mrna (follows demo)
#' @param dataset : dataset that needs to consolidated (follows format demo)
#' @param PC_VAR_THRESH (0.09): The % threshold for PCs to accept based on scale of eigenvalue (must be <1)
#' @param pc_collapse : The additional consolidation function applied to selected PCs. (default take the max)
#' @param DEBUG FALSE: flag to print debug statements
#'
#' @export
dataset.collapse.pc.singular <- function(mrna,
                                         dataset,
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

  # Note: parallel later ...
  data.collapsed <- sapply(1:p,FUN=function(i){collapse_gene_data(i,dataset)})
  row.names(data.collapsed) <- row.names(mrna)
  colnames(data.collapsed) <- colnames(mrna)
  return(data.collapsed)
}

#' get.data.names
#' @description get a vector of usable names for the list of datasets.
get.data.names <- function(data.size,
                           data.names,
                           sep="_",
                           default.data.name="data",
                           DEBUG = FALSE){
  if(DEBUG){
    print("in get.data.names")
    print(data.names)
  }
  if(is.null(data.names)){
    data.names <- rep(default.data.name, data.size)
  }
  data.names <- unname(sapply(data.names, FUN = function(name){ifelse(name == "",default.data.name,paste(name))}))
  if(DEBUG){
    print(data.names)
  }
  name.counter <- list()
  for(i in 1:data.size){
    if(DEBUG){
      print("looping over data.names")
      print(data.names)
      print(i)
      print(name.counter)
    }
    if(sum(data.names[i] == data.names) > 1 || exists(data.names[i],name.counter)){
      if(exists(data.names[i],name.counter)){
        name.counter[data.names[i]] <- as.integer(name.counter[data.names[i]]) + 1
      } else {
        name.counter[data.names[i]] <- 1
      }
      new.name <- sprintf("%s%s%d",data.names[i],sep,as.integer(name.counter[data.names[i]]))
      if(new.name %in% data.names[-1*i]){
        stop(paste("collision on data renaming '", new.name, "' Check sep.", sep="", collapse = ""))
      }
      data.names[i] <- new.name
    }
  }
  return(data.names)
}
