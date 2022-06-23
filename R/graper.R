#' @title Fit a regression model with graper
#' @name graper
#' @description Fit a regression model with graper given a matrix
#'  of predictors (\code{X}), a response vector (\code{y}) and
#'  a vector of group memberships for each predictor
#'  in \code{X} (\code{annot}).
#'  For each group a different strength of penalization
#'  is determined adaptively.
#' @param X design matrix of size n (samples) x p (features)
#' @param y response vector of size n
#' @param annot factor of length p indicating group membership
#'  of each feature (column) in X
#' @param family Likelihood model for the response,
#'  either "gaussian" for linear regression or
#'  "binomial" for logistic regression
#' @param factoriseQ if set to TRUE, the variational
#'  distribution is assumed to fully factorize across
#'  features (faster, default). If FALSE, a
#'  multivariate variational distribution is used.
#' @param spikeslab if set to TRUE, a spike and slab prior
#'  on the coefficients (default).
#' @param d_tau hyper-parameters for prior of tau (noise precision)
#' @param r_tau hyper-parameters for prior of tau (noise precision)
#' @param d_gamma hyper-parameters for prior of gamma
#'  (coefficients' prior precision)
#' @param r_gamma hyper-parameters for prior of gamma
#'  (coefficients' prior precision)
#' @param r_pi hyper-parameters for Beta prior of the mixture
#'  probabilities in the spike and slab prior
#' @param d_pi hyper-parameters for Beta prior of the mixture
#'  probabilities in the spike and slab prior
#' @param max_iter maximum number of iterations
#' @param th convergence threshold for the evidence lower bound (ELB)
#' @param intercept whether to include an intercept into the model
#' @param calcELB whether to calculate the evidence lower bound (ELB)
#' @param verbose  whether to print out intermediate messages
#' during fitting
#' @param freqELB frequency at which the evidence lower bound (ELB)
#' is to be calculated, i.e. each freqELB-th iteration
#' @param n_rep number of repetitions with different random
#' initializations  to be fit
#' @param standardize whether to standardize the predictors
#' to unit variance
#' @param init_psi initial value for the spike variables
#' @param nogamma if TRUE, the normal prior will have same
#'  variance for all groups (only relevant for spikeslab = TRUE)
#' @details The function trains the graper model given
#' a matrix of predictors (\code{X}), a response vector (\code{y})
#' and a vector of group memberships for each predictor in \code{X}
#' (\code{annot}). For each feature group as specified in
#' \code{annot} a penalty factor and sparsity level is learnt.
#'
#' By default it uses a Spike-and-Slab prior on the coefficients
#' and uses a fully factorized variational distribution
#' in the inference. This provides a fast way to train the model.
#' Using \code{spikeslab=FALSE} a ridge regression like model can
#' be fitted using a normal instead of the spike and slab prior.
#' Setting \code{factoriseQ = FALSE} gives a more exact inference
#' scheme based on a multivariate variational distribution,
#' but can be much slower.
#'
#'  As the optimization is non-convex is can be helpful to
#'  use multiple random initializations by setting \code{n_rep}
#'  to a value larger 1. The returned model is then chosen
#'  as the optimal fit with respect to the evidence lower bound (ELB).
#'
#'  Depending on the response vector a linear regression model
#'  (\code{family = "gaussian"}) or a logistic regression model
#'  (\code{family = "binomial"}) is fitted.
#'  Note, that the implementation of logistic regression is still
#'  experimental.
#'
#' @return A graper object containing
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
#' @import Rcpp
#' @importFrom matrixStats colSds
#' @export
#' @examples
#' # create data
#' dat <- makeExampleData()
#'
#' # fit a sparse model with spike and slab prior
#' fit <- graper(dat$X, dat$y, dat$annot)
#' fit # print fitted object
#' beta <- coef(fit, include_intercept=FALSE) # model coeffients
#' pips <- getPIPs(fit) # posterior inclusion probabilities
#' pf <- fit$EW_gamma # penalty factors per group
#' sparsities <- fit$EW_pi # sparsity levels per group
#'
#' # fit a dense model without spike and slab prior
#' fit <- graper(dat$X, dat$y, dat$annot, spikeslab=FALSE)
#'
#' # fit a dense model using a multivariate variational distribution
#' fit <- graper(dat$X, dat$y, dat$annot, factoriseQ=TRUE,
#'       spikeslab=FALSE)


graper <- function(X, y, annot, factoriseQ = TRUE, spikeslab = TRUE,
    intercept = TRUE, family = "gaussian", standardize = TRUE, n_rep = 1,
    max_iter = 3000, th = 0.01, d_tau = 0.001, r_tau = 0.001,
    d_gamma = 0.001, r_gamma = 0.001, r_pi = 1, d_pi = 1, calcELB = TRUE,
    verbose = TRUE, freqELB = 1, nogamma = FALSE, init_psi = 1) {

  stopifnot(ncol(X) == length(annot)) # check correct dimensions of input
  if(!spikeslab & !nogamma) nogamma <- FALSE # nogamma only with spikeslab
  annot <- factor(annot, levels=unique(annot))
  p <- ncol(X)  # no. of features
  g <- length(unique(annot)) # no. of groups
  NoPerGroup <- vapply(unique(annot), function(x){ # features per group
    sum(annot == x)}, numeric(1))
  names(NoPerGroup) <- unique(annot)
  if(verbose){
    message("Fitting a model with ", g, " groups, ", nrow(X),
      " samples and ", p , " features.")
  }
  
  calcELB <- .checkCalcELB(calcELB=calcELB, family=family)
  n_rep <- .checkNrep(n_rep=n_rep, calcELB=calcELB)
  if(standardize) {
    sf <-  matrixStats::colSds(X) # scale by sd
    X <- scale(X, center = FALSE, scale=sf)
  } else sf <- rep(1,p)

  reslist <- lapply(seq_len(n_rep), function(rep){
    .graperSingle(rep, family, X, y, annot, factoriseQ, spikeslab,
      intercept, max_iter, th, d_tau, r_tau, d_gamma, r_gamma, r_pi, d_pi,
      calcELB, verbose, freqELB, nogamma, init_psi, p, g, NoPerGroup)
  })
  res <- .selectBestModel(reslist=reslist, n_rep=n_rep) # select best model
  if(standardize) { # revert coefficients to original scale
    res$EW_beta <- res$EW_beta / sf
    if(!factoriseQ) {
      res$Sigma_beta <- diag(1 / sf) %*% res$Sigma_beta %*% diag(1 / sf)
    } else {
      res$Sigma_beta <- diag(1/(sf^2) * diag(as.matrix(res$Sigma_beta)))
    }
  }
  res$annot <- annot
  res$Options <- list(factoriseQ = factoriseQ, spikeslab = spikeslab,
    d_tau = d_tau, r_tau = r_tau, d_gamma = d_gamma, r_gamma = r_gamma,
    r_pi = r_pi, d_pi = d_pi, max_iter = max_iter, th = th,
    intercept = intercept, calcELB = calcELB, verbose = verbose,
    freqELB = freqELB, family = family, nogamma = nogamma,
    standardize = standardize, featurenames = colnames(X))
  if(all(is.na(res$ELB_trace) | !is.finite(res$ELB_trace))) {
    res$ELB_trace <- NULL # remove ELB slot if not calculated
  }
  class(res) <- "graper"
  return(res)
}

# function to fit the graper model for a Gaussian response
.graperGaussian <- function(X, y, annot, factoriseQ, spikeslab,
        intercept, max_iter, th, d_tau, r_tau, d_gamma,
        r_gamma, r_pi, d_pi, calcELB, verbose, freqELB,
        nogamma, init_psi, p, g, NoPerGroup) {
    if (intercept) { # remove intercept effect by centering X and y
      X <- scale(X, center = TRUE, scale = FALSE)
      y <- scale(y, center = TRUE, scale = FALSE)
    }
    if (spikeslab) {
      if (!factoriseQ) {
        factoriseQ <- TRUE
        warning("Seeting factoriseQ to TRUE")
      }
      # initialize slab mean and spike prob.
      mu_init <- rnorm(p)
      psi_init <- rep(init_psi,p)
      if(!nogamma) {
        res <- graperCpp_sparse_ff(X, y, annot, g, NoPerGroup, d_tau,
            r_tau, d_gamma, r_gamma, r_pi, d_pi, max_iter, th,
            calcELB, verbose, freqELB, mu_init, psi_init)
      } else {
        res <- graperCpp_sparse_ff_nogamma(X, y, annot, g,
          NoPerGroup, d_tau, r_tau, d_gamma, r_gamma, r_pi, d_pi,
          max_iter, th, calcELB, verbose, freqELB, mu_init, psi_init)
      }
    } else {
      if (factoriseQ) {
        mu_init <- rnorm(p) # initialize randomly
        res <- graperCpp_dense_ff(X, y, annot, g, NoPerGroup,
            d_tau, r_tau, d_gamma, r_gamma, max_iter, th,
            calcELB, verbose,freqELB, mu_init)
      } else {
        message("factoriseQ = FALSE might take some time
          to compute. Set factoriseQ = TRUE for fast solution.")
        res <- graperCpp_dense_nf(X, y, annot, g, NoPerGroup, d_tau,
            r_tau, d_gamma, r_gamma, max_iter, th, calcELB,
            verbose, freqELB)
      }
    }
    res$intercept <- NULL
    if (intercept) { # calculate intercept
      res$intercept <- attr(y, "scaled:center") -
          sum(attr(X, "scaled:center") * res$EW_beta)
    }
    if(!nogamma) {
      rownames(res$EW_gamma) <- unique(annot)
    }
    res
}

# function to fit the graper model for a Bernoulli response
.graperBinomial <- function(X, y, annot, factoriseQ, spikeslab,
        intercept, max_iter, th, d_gamma, r_gamma,
        r_pi, d_pi, calcELB, verbose, freqELB,
        init_psi, p, g, NoPerGroup) {
      if (spikeslab){
        if (!factoriseQ) {
          factoriseQ <- TRUE
          warning("Using fully factorized approach
              with a spike and slab prior")
        }
        # initialize slab mean and spike prob. randomly
        mu_init <- rnorm(p)
        # psi_init <- runif(p)
        psi_init <- rep(init_psi,p)
        res <- graperCpp_sparse_logistic_ff(X, y, annot,
            g, NoPerGroup, d_gamma, r_gamma,
            r_pi, d_pi, max_iter, th, calcELB,verbose,
            freqELB, mu_init, psi_init, intercept)
      } else {
        if (factoriseQ){
          # initialize coefficients mean randomly
          mu_init <- rnorm(p)
          res <- graperCpp_logistic_ff(X, y, annot,
              g, NoPerGroup, d_gamma, r_gamma,
              max_iter, th, calcELB, verbose, freqELB,
              mu_init, intercept)
        } else {
          warning("factoriseQ=FALSE is not maintained currently
              for the logistic model.
              No intercept option and ELB available.")
          res <- graperCpp_logistic_nf(X, y, annot,
              g, NoPerGroup, d_gamma, r_gamma,
              max_iter, th, calcELB, verbose, freqELB)
        }
      }

      if(!intercept) {
        res$intercept <- NULL
      }

      res
}

# function to select the best model (by the ELBO)
.selectBestModel <- function(reslist, n_rep) {
  if(n_rep == 1) {
    res <- reslist[[1]]
  } else {
    best_idx <- which.max(vapply(reslist, function(l) l$ELB, numeric(1)))
    if(is.na(best_idx) | is.null(best_idx)) {
      warning("Model selection based on ELB encountered errors.
          Returned model is picked arbitrarily!")
      best_idx <- 1
    }

    res <- reslist[[best_idx]]
  }
}

# function to fit a single graper model
.graperSingle <- function(rep, family, X, y, annot, factoriseQ, spikeslab,
        intercept, max_iter, th, d_tau, r_tau, d_gamma,
        r_gamma, r_pi, d_pi, calcELB, verbose, freqELB,
        nogamma, init_psi, p, g, NoPerGroup) {
    if(verbose){
      message("Fitting with random init ", rep)
    }

    if (family == "gaussian") {
      .graperGaussian(X=X, y=y, annot=annot,
        factoriseQ=factoriseQ, spikeslab=spikeslab,
        intercept=intercept, max_iter=max_iter, th=th,
        d_tau=d_tau, r_tau=r_tau, d_gamma=d_gamma,
        r_gamma=r_gamma, r_pi=r_pi, d_pi=d_pi,
        calcELB=calcELB, verbose=verbose, freqELB=freqELB,
        nogamma=nogamma, init_psi=init_psi, p=p, g=g,
        NoPerGroup=NoPerGroup)
    } else if (family == "binomial") {
      .graperBinomial(X=X, y=y, annot=annot,
        factoriseQ=factoriseQ, spikeslab=spikeslab,
        intercept=intercept, max_iter=max_iter, th=th,
        d_gamma=d_gamma, r_gamma=r_gamma, r_pi=r_pi,
        d_pi=d_pi, calcELB=calcELB, verbose=verbose,
        freqELB=freqELB, init_psi=init_psi,
        p=p, g=g, NoPerGroup=NoPerGroup)
    } else {
      stop("Family not implemented.
        Needs to be either binomial or gaussian.")
    }
}

# function to check ELB argument is ok
.checkCalcELB <- function(calcELB, family){
  if(family == "binomial" & calcELB){
    calcELB <- FALSE
    warning("The implementation of logistic regression
      is still experimental, ELB is not calculated.")
  }
  calcELB
}

# function to check n_rep argument is ok
.checkNrep <- function(n_rep, calcELB){
  if(!calcELB & n_rep >1) {
    warning("Only using a single trial now as calcELB = FALSE.")
    n_rep <- 1
  }
  n_rep
}
