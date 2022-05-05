# objects.R
# this file consists of the main iBAG object definition and manipulators.
# this also includes all other supplementary objects that the iBAG obj uses.

# required to build the classes
library(R6)

#' iBAG_data
#'
#' @name iBAG_data
#' @description R6 class to store and manage the iBAG data.
#' This does not contain or manage model result data.
#'
iBAG_data <- R6Class("iBAG_data",
  public = list(
    #' @description constructor
    #' @usage iBAG_data$new()
    #'
    #' @param mrna (iBAG::demo_mrna) dataframe of mrna data
    #' @param outcome (iBAG::demo_outcome) dataframe of outcome data
    #' @param data (list(cnv = iBAG::demo_cnv)) list of upstream data
    #' @param DEBUG (FALSE) initialize object in DEBUG mode
    #' @param validate (TRUE) to validate the data supplied. Will raise Error if data is not valid.
    #' @param one_val_per_gene (TRUE) whether all the upstream data has 1 column per gene or not.
    #' @param ... : ...
    initialize = function(mrna = iBAG::demo_mrna,
                          outcome = iBAG::demo_outcome,
                          data = list(cnv=iBAG::demo_cnv,meth=iBAG::demo_meth),
                          DEBUG = FALSE,
                          validate = TRUE,
                          one_val_per_gene = TRUE,...){
      private$DEBUG = DEBUG
      private$one_val_per_gene = one_val_per_gene
      kwargs <- list(...)
      if(validate){
        sep = "_"
        if(one_val_per_gene){
          sep = NULL
        } else if(!is.null(kwargs$sep)){
          sep = kwargs$sep
        }
        if(!iBAG::data.validate(mrna,outcome,data,sep,DEBUG)){
          stop("Data is not in valid format. Please consult ?iBAG::data.validate or demo_data for examples of valid data")
        }
      }
      private$mrna = mrna
      private$outcome = outcome
      private$data = data

      # set inferred data
      private$genes <- colnames(mrna)
      private$patients <- row.names(mrna)
      private$data_names <- names(data)
      private$n_genes <- length(private$genes)
      private$n_patients <- length(private$patients)
      private$n_data <- length(private$data_names)
    },
    #' @description get the dataset specified by user
    #' @usage iBAG_data$new()$get.data()
    #'
    #' @param index (NULL): index of the dataset to return. Can either be an integer otherwise it's treated as a string. Can be null (sends entire list back).
    #' @return dataframe of the requested dataset or list of all dataframes
    #'
    #' TODO: address not found/out of bounds conditions
    get.data = function(index = NULL){
      if(private$DEBUG){
        print(index)
        print(dim(private$data[[index]]))
      }
      if(is.null(index)){
        return(private$data)
      } else if(is.numeric(index) && index == as.integer(index)){
        return(private$data[[index]])
      } else {
        index = paste(index)
        return(private$data[index][[1]])
      }
    },
    #' @description get the mrna dataset
    #' @usage iBAG_data$new()$get.mrna()
    #'
    #' @return dataframe
    get.mrna = function(){
      if(private$DEBUG){
        print(dim(private$mrna))
      }
      return(private$mrna)
    },
    #' get.outcome
    #'
    #' @description get the outcome dataset
    #' @usage iBAG_data$new()$get.outcome()
    #'
    #' @return dataframe
    get.outcome = function(){
      if(private$DEBUG){
        print(dim(private$outcome))
      }
      return(private$outcome)
    },
    #' @description get the number of upstream datasets
    #' @usage iBAG_data$new()$get.n_data()
    #'
    #' @return numeric (integer)
    get.n_data = function(){
      return(length(private$data))
    },
    #' get.n_patients
    #'
    #' @description get the number of patients in dataset
    #' @usage iBAG_data$new()$get.n_patients()
    #'
    #' @return numeric (integer)
    get.n_patients = function(){
      return(dim(private$mrna)[1])
    },
    #' @description get a vector of patients in dataset
    #' @usage iBAG_data$new()$get.patients()
    #'
    #' @return vector of numeric or strings
    get.patients = function(){
      return(row.names(private$mrna))
    },
    #' get.n_genes
    #'
    #' @description get the number of genes in dataset
    #' @usage iBAG_data$new()$get.n_genes()
    #'
    #' @return numeric (integer)
    get.n_genes = function(){
      return(dim(private$mrna)[2])
    },
    #' @description get a vector of genes in dataset
    #' @usage iBAG_data$new()$get.genes()
    #'
    #' @return vector of numeric or strings
    get.genes = function(){
      return(colnames(private$mrna))
    }
  ),
  private = list(
    # roxygen was being pissy about documenting private variables, but I left em in ;)
    # @field data list of dataframes for the upstream data
    data = NULL,
    # @field mrna dataframe for the mrna data
    mrna = NULL,
    # @field outcome dataframe for the outcome data
    outcome = NULL,
    # @field DEBUG boolean flag used to figure out what went wrong or print stuff to console
    DEBUG = NULL,
    # @field one_val_per_gene boolean flag used to define whether there's only 1 column of data per gene in upstream data or not (used to speed up computations).
    one_val_per_gene = NULL,
    # @field genes vector of gene names (taken from colnames of mrna)
    genes = NULL,
    # @field patients vector of patient names (taken from row.names of mrna)
    patients = NULL,
    # @field data_names vector of data source names (taken from names(data))
    data_names = NULL,
    # @field n_genes number of genes (length(genes))
    n_genes = NULL,
    # @field n_patients number of patients (length(patients))
    n_patients = NULL,
    # @field n_data number of data (length(data_names))
    n_data = NULL
  )
)

#' iBAG_results
#'
#' @name iBAG_results
#' @description R6 class to store and manage the iBAG results.
#' This does not contain a copy of the iBAG_data object.
#'
iBAG_results <- R6Class("iBAG_results",public = list(),private = list())
