# objects.R
# this file consists of the main iBAG object definition and manipulators.
# this also includes all other supplementary objects that the iBAG obj uses.

# required to build the classes
library(R6)

#' iBAG_data
#'
#' @name iBAG_data
#'
#' @export
iBAG_data <- R6Class("iBAG_data",
  public = list(
    #'
    #' @param DEBUG boolean: initialize object in DEBUG mode
    #' @param print_status boolean: initialize object to print (to stdout) everything that it does
    #' @param ... : ...
    #'
    initialize = function(DEBUG = FALSE,print_status=FALSE,...){
      # This constructor is quite flexible
      # meth_file character: file name for the methylation data file
      # mrna_file character: file name for the mrna expression data file
      # cnv_file character: file name for the cnv data file
      # outcome_file character: file name for the outcome file
      # meth dataframe: methylation data (takes precedence over meth_file)
      # mrna dataframe: mrna data (takes precedence over the mrna_file)
      # cnv dataframe: cnv data (takes precedence over cnv_file)
      # outcome dataframe: outcome data (takes precedence over the outcome_file)
      # DEBUG boolean = FALSE: flag for whether to print DEBUG information
      # print_status boolean = FALSE: flag for whether to print progress

      # status flags
      private$DEBUG = DEBUG
      private$print_status = print_status

      #
      #private$meth_file = system.file("data","methylationdata.csv",package = "iBAGpkg")
      #private$mrna_file = system.file("data","mrnadata.csv",package = "iBAGpkg")
      #private$cnv_file = system.file("data","copynumberdata.csv",package = "iBAGpkg")
      #private$outcome_file = system.file("data","survivaltimes.csv",package = "iBAGpkg")

      private$meth = iBAGpkg::demo_meth
      private$mrna = iBAGpkg::demo_mrna
      private$cnv = iBAGpkg::demo_cnv
      private$outcome = iBAGpkg::demo_outcome

      kwargs <- list(...)
    },
    #'
    get.meth = function(){
      return(private$meth)
    },
    #'
    get.cnv = function(){
      return(private$cnv)
    },
    #'
    get.mrna = function(){
      return(private$mrna)
    },
    #'
    get.outcome = function(){
      return(private$outcome)
    }
  ),
  private = list(
    n_genes = NA,
    n_patients = NA,
    mrna = NULL,
    meth = NULL,
    cnv = NULL,
    outcome = NULL,
    DEBUG = NULL,
    print_status = NULL
  )
)
