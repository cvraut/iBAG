# clin.R
# This file contains the functions to fit the clinical model

#' clin.model
#'
#' @description This fits the clinical model for iBAG
#' 
#' @details The clinical model is based on a spike-and-slab prior for Bayesian Linear Regression.
#' The function uses Variational Bayes by default to fit the model. Implementation follows the 
#' graper package by Brita Valten & is embedded & slightly modified in this package.
#' 
#' The function primarily takes in data from the output of mech.model and construct.priors.
#' 
#' @param X the covariate matrix (use result from mech.model)
#' @param Y the vector of output (use result from mech.model)
#' @param priors (iBAG::construct.priors()): the list of priors , uses non-informative priors by default
#' @param DEBUG (FALSE): debug flag
#' @param ... extra arguments for the iBAG model (passed to graper function)
#' @return graper object containing
#' \describe{
#' \item{EW_beta}{estimated model coefficients in liner/logistic
#' regression}
#' \item{EW_s}{estimated posterior-inclusion probabilities
#'  for each feature}
#' \item{intercept}{estimated intercept term}
#' \item{annot}{annotation vector of features to the groups as
#'  specified when calling \code{\link{graper}}}
#' \item{EW_gamma}{estimated penalty factor per group}
#' \item{EW_pi}{estimated sparsity level per group
#'  (from 1 (dense) to 0 (sparse))}
#' \item{EW_tau}{estimated noise precision}
#' \item{sigma2_tildebeta_s1, EW_tildebeta_s1, alpha_gamma,
#'  alpha_tau, beta_tau, Sigma_beta, alpha_pi, beta_pi}{parameters
#'  of the variational distributions of beta, gamma, tau and pi}
#' \item{ELB}{final value of the evidence lower bound}
#' \item{ELB_trace}{values of the  evidence lower bound
#'  for all iterations}
#' \item{Options}{other options used when calling \code{\link{graper}}}
#' }
#' @export
clin.model <- function(X,
                       Y,
                       priors = iBAG::construct.priors(),
                       DEBUG = FALSE,
                       ...){
  graper.args <- c(list(X = X, y = Y, annot = rep(1, ncol(X))), priors, list(...))

  # return(iBAG::graper(X = X, 
  #                     y = Y, 
  #                     annot = rep(1,ncol(X)), 
  #                     d_tau = priors$d_tau,
  #                     r_tau = priors$r_tau,
  #                     d_pi = priors$d_pi,
  #                     r_pi = priors$r_pi,
  #                     r_gamma = priors$r_gamma,
  #                     d_gamma = priors$d_gamma,
  #                     ...)))
  return(do.call(iBAG::graper,graper.args))
}