# priors.R
# file containing prior parameter setting for the clinical model


#' construct.priors
#' @name construct.priors
#' @description constructs the priors for the clinical model. includes basic bound checks
#' @details
#' follows from https://arxiv.org/pdf/1811.02962v1.pdf section 2.3
#' - (residual precision) tau ~ Gamma(r_tau,d_tau)
#' - (prob of kth coefficient != 0) pi_k ~ Beta(d_pi, r_pi)
#' - (precision of kth coefficient) gamma_k ~ Gamma(r_gamma,d_gamma)
#' Note: all parameters must be positive (>0)
#' @usage construct.priors()
#' 
#' @param r_tau (0.001) hyperparameter for residual precision tau
#' @param d_tau (0.001) hyperparameter for residual precision tau
#' @param d_pi (1) hyperparameter for prob of kth coefficient != 0
#' @param r_pi (1) hyperparameter for prob of kth coefficient != 0
#' @param r_gamma (0.001) hyperparameter for precision of kth coefficient
#' @param d_gamma (0.001) hyperparameter for precision of kth coefficient
#' @return list with the following hyperparamenters accessible by index or param name
#' 
#' 
construct.priors <- function(r_tau = 0.001,
                             d_tau = 0.001,
                             d_pi = 1,
                             r_pi = 1,
                             r_gamma = 0.001,
                             d_gamma = 0.001){
  priors <- list(r_tau = r_tau,
                 d_tau = d_tau,
                 d_pi = d_pi,
                 r_pi = r_pi,
                 r_gamma = r_gamma,
                 d_gamma = d_gamma)
  if(any(priors <= 0)){
    stop(sprintf("All prior params must be > 0. supplied params: %s", toString(priors)))
  }
  return(priors)
}