library(iBAG)

testthat::test_that("test graper works", {
  set.seed(19890421)
  # build graper model from following params
  n <- 100
  pg <- c(100, 100, 10, 10)
  gammas <- c(0.1, 10, 0.1, 10)
  pis <- c(0.5, 0.5, 0.5, 0.5)
  tau <- 1
  rho <- 0.1
  response <- "gaussian"
  intercept <- 0

  g <- length(pg)
  p <- sum(pg)

  # construct design
    Sigma <- stats::toeplitz(rho ^ (0:(p - 1)))
    X <- matrix(stats::rnorm(n * p), n, p) %*% chol(Sigma)
    X <- scale(X)

    # simulate coefficients
    beta_tilde <- Reduce(c, lapply(seq_len(g), function(k) {
        stats::rnorm(pg[k], 0, sqrt(1 / gammas[k]))
    }))
    s <- Reduce(c,lapply(seq_len(g), function(k) {
        stats::rbinom(pg[k], 1, pis[k])
    }))
    beta <- s * beta_tilde

    #simulate response
    if(response == "gaussian"){
        y <- stats::rnorm(n, X %*% beta + intercept, 1 / sqrt(tau))
    } else {
        y <- stats::rbinom(n, 1, 1 / (1 + exp(- (X %*% beta + intercept))))
    }

    #group annotations
    annot <- rep(seq_along(pg), times=pg)

    graper.result <- iBAG::graper(X,y,annot,max_iter=10000,n_rep=5,verbose = FALSE)

    #testthat::expect_equal(beta,graper.result$EW_beta,tolerance = 0.01)
    testthat::expect_equal(c(beta), c(graper.result$EW_beta), tolerance = 2,
    label = "test that the estimates for beta are reasonably close within 2 of the true")
})