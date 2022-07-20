# iBAG.R
# This file contains roxygen2 documentation for the pkg. Used primarily for NAMESPACE construction
# This file also contains the main iBAG() function that brings it all together

#' iBAG: A package for fitting the iBAG model.
#'
#' Get started with iBAG in 2 lines:
#' demo_data <- iBAG::demo_data()
#' iBAG_result <- iBAG::iBAG(demo_data)
#' summary(iBAG_result)
#' summary(iBAG_result$mechmodel)
#' 
#' @section iBAG functions:
#' The iBAG functions ...
#'
#' @docType package
#' @name iBAG
#' @useDynLib iBAG
NULL
#> NULL

#' iBAG
#' 
#' @name iBAG
#' @description fits the iBAG model
#' @details 
#' follows details from ____
#' 
#' consult [supplementary] for an extended explanation & example
#' @usage 
#' data <- iBAG::demo_data()
#' result <- iBAG::iBAG(data)
#' summary(result)
#' 
#' @param ... additional kwargs
#' @return iBAG object -> main components are iBAG_object$mechmodel and iBAG_object$clinmodel
#' @export 
iBAG <- function(...){
  return()
}