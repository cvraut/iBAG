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
    initialize = function(mrna = iBAG::demo_mrna,
                          outcome = iBAG::demo_outcome,
                          data1 = iBAG::demo_cnv,
                          DEBUG = FALSE,
                          print_status=FALSE,...){
      private$DEBUG = DEBUG
      private$print_status = print_status

      private$mrna = mrna
      private$outcome = outcome
      private$data = list(data1)

      kwargs <- list(...)
    },
    #'
    get.data = function(k=1){
      if(private$DEBUG){
        print(k)
        print(dim(private$data[[k]]))
      }
      return(private$data[[k]])
    },
    #'
    get.mrna = function(){
      if(private$DEBUG){
        print(dim(private$mrna))
      }
      return(private$mrna)
    },
    #'
    get.outcome = function(){
      if(private$DEBUG){
        print(dim(private$outcome))
      }
      return(private$outcome)
    },
    get.k = function(){
      return(length(private$data))
    },
    n_patients = function(){
      return(dim(private$mrna)[1])
    },
    get.patients = function(){
      return(row.names(private$mrna))
    },
    n_genes = function(){
      return(dim(private$mrna)[2])
    },
    get.genes = function(){
      return(colnames(private$mrna))
    }
  ),
  private = list(
    data = NULL,
    mrna = NULL,
    outcome = NULL,
    DEBUG = NULL,
    print_status = NULL
  )
)
