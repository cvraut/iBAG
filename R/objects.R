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
                          data = list(cnv=iBAG::demo_cnv, meth=iBAG::demo_meth),
                          DEBUG = FALSE,
                          validate = TRUE,
                          one_val_per_gene = TRUE, ...){
      private$DEBUG <- DEBUG
      private$one_val_per_gene <- one_val_per_gene
      kwargs <- list(...)
      if(validate){
        sep <- "_"
        if(one_val_per_gene){
          sep <- NULL
        } else if(!is.null(kwargs$sep)){
          sep <- kwargs$sep
        }
        if(!iBAG::data.validate(mrna, outcome, data, sep, DEBUG)){
          stop(paste("Data is not in valid format.",
          "Please consult ?iBAG::data.validate or demo_data for examples of valid data"), sep="\n")
        }
      }
      private$mrna <- mrna
      private$outcome <- outcome
      private$data <- data

      # set inferred data
      private$genes <- colnames(mrna)
      private$patients <- row.names(mrna)
      private$data_names <- names(data)
      private$n_genes <- length(private$genes)
      private$n_patients <- length(private$patients)
      private$n_data <- length(private$data_names)
    },
    #' get.data
    #' 
    #' @description get the dataset specified by user
    #' @usage iBAG_data$new()$get.data()
    #'
    #' @param index (NULL): index of the dataset to return. Can either be integer or string. if null, return list.
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
        index <- paste(index)
        return(private$data[index][[1]])
      }
    },
    #' get.mrna
    #' 
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
    #' get.n_data
    #' 
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
    #' get.genes
    #' 
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
    # @field one_val_per_gene boolean flag whether there's only 1 column of data per gene in upstream
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
iBAG_results <- R6Class("iBAG_results",
  public = list(
    #' constructor
    #' @description constructor
    #' @usage iBAG_results$new()
    #'
    #' @param X (NULL): the total X matrix generated by the mech model
    #' @param Y (NULL): the vector of outcomes that goes into the clinical model
    #' @param SS (list(NULL,NULL)): a list of vectors representing the Sum of Squares Error from the mechanistic model.
    #' @param beta_mean (NULL): a vector of posterior means for the coefficients from the clinical model.
    #' @param beta_incl_prob (NULL): a vector of posterior inclusion probabilities for the beta coefficients.
    #' @param DEBUG (FALSE): initialize object in DEBUG model.
    #' @param validate (TRUE): validate the data coming in from the constructor (or any of the setters)
    #' @param ... : ...
    #'
    #' @details This constructs a class containing all the results from the iBAG object.
    #' Some additional details:
    #'  - SS:
    #'    - The list size is always >2. The 1st element is the SST, the last element is the SSO,
    #'    the middle elements correspond to the sum of squares from the datasets.
    #'  - validate:
    #'    - minimum data required to validate is X,Y,SS (all produced by mechmodel)
    #'      - if this data is not provided will issue warning & set validate to FALSE
    #'    - n_data is derived from length(SS)-2
    #'    - n_patients is derived as the nrow(X)
    #'    - check that ncol(X)%n_data == 0
    #'      - if not throw error
    #'      - if true: n_genes = ncol(X)%/%n_data
    #'    - check that length(Y) == n_patients
    #'      - if not throw error
    #'    - for each item in SS
    #'      - check that length(item) == n_genes
    #'        - if not throw error
    #'    - if beta_mean & beta_incl_prob are included
    #'      - check that their length == ncol(X)
    #'        - if not throw error
    #'    - This does not check that the data is in a numeric form! (Check that yourself)
    initialize = function(X = matrix(0, 1, 1),
                          Y = c(0),
                          SS = list(c(0), c(0), c(0)),
                          beta_mean = c(0),
                          beta_incl_prob = c(0),
                          DEBUG = FALSE,
                          validate = TRUE,
                          ...){
      private$DEBUG <- DEBUG
      private$validate <- validate
      
      private$n_data <- length(SS) - 2
      private$n_patients <- nrow(X)
      private$n_genes <- ncol(X)%/%private$n_data
      
      if(is.null(X) || is.null(Y) || is.null(SS)){
        warning("X, Y, or SS was missing. Setting validate to FALSE")
        validate <- FALSE
      } else if(validate){
        if(length(Y) != private$n_patients){
          stop("n_patients is not consistent across X and Y")
        } else if(private$n_data <= 0){
          stop("SS does not contain at least 3 elements")
        } else if(ncol(X) != private$n_genes * private$n_data){
          stop("X does not contain the right amount of gene columns according to SS")
        } else if(any(sapply(SS, FUN=length) != private$n_genes)){
          stop("SS does not contain the correct amount of data for each gene")
        } else if(!is.null(beta_mean) && length(beta_mean) != ncol(X)){
          stop("dimensions of beta_mean not consistent with X")
        } else if(!is.null(beta_incl_prob) && length(beta_incl_prob) != ncol(X)){
          stop("dimensions of beta_incl_prob not consistent with X")
        }
      }
      private$X <- X
      private$Y <- Y
      private$SS <- SS
      private$beta_mean <- beta_mean
      private$beta_incl_prob <- beta_incl_prob
      private$validate <- validate
    },
    #' get.X
    #' 
    #' @description returns X matrix
    #' 
    #' @return matrix that should be n_patients rows and n_genes*n_data columns
    get.X = function(){
      return(private$X)
    },
    #' get.Y
    #' 
    #' @description returns Y vector
    #' 
    #' @return vector that is n_patients long
    get.Y = function(){
      return(private$Y)
    },
    #' get.SS
    #' 
    #' @description returns the Sum of Squares mech model results, either all or some of them
    #' 
    #' @param index (NULL) can be NULL, integer, or string. If NULL returns list, otherwise returns the vector at index
    #' @return a vector or a list of vectors 
    get.SS = function(index = NULL){
      if(is.null(index)){
        return(private$SS)
      } else if(is.numeric(index) && index == as.integer(index)){
        return(private$SS[[index]])
      } else {
        index <- paste(index)
        return(private$SS[index][[1]])
      }
    },
    #' get.beta_mean
    #' 
    #' @description returns the vector of posterior beta means
    #' 
    #' @return a vector of posterior beta means
    get.beta_mean = function(){
      return(private$beta_mean)
    },
    #' get.beta_incl_prob
    #' 
    #' @description returns a vector of posterior beta inclusion probabilities
    #' 
    #' @return a vector of posterior beta inclusion probabilities
    get.beta_incl_prob = function(){
      return(private$beta_incl_prob)
    },
    #' get.n_patients
    #' 
    #' @description returns the number of patients
    #' 
    #' @return integer of number of patients
    get.n_patients = function(){
      return(private$n_patients)
    },
    #' get.genes
    #' 
    #' @description returns the number of genes
    #' 
    #' @return integer of number of genes
    get.n_genes = function(){
      return(private$n_genes)
    },
    #' get.n_data
    #' 
    #' @description returns the number of upstream datasets
    #' 
    #' @return integer of number of upstream datasets
    get.n_data = function(){
      return(private$n_data)
    }
  ),
  private = list(
    X = NULL,
    Y = NULL,
    SS = NULL,
    beta_mean = NULL,
    beta_incl_prob = NULL,
    DEBUG = NULL,
    validate = NULL,
    n_patients = NULL,
    n_genes = NULL,
    n_data = NULL
  )
)