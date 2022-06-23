library(iBAG)

testthat::test_that("test clin.model",{
  set.seed(19890421)
  # test on dummy data
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

    graper.result <- iBAG::clin.model(X,y,max_iter=50,verbose = FALSE)
    testthat::expect_equal(graper.result$ELB,-487,tolerance = 1,label = "checking clin.model via ELB on dummy data")

    # test on the demo data
    demo_data <- iBAG::iBAG_data$new()
    mech.res <- iBAG::mech.model(demo_data$get.mrna(),demo_data$get.data())
    testthat::expect_silent(clin.res <- iBAG::clin.model(mech.res$X,demo_data$get.outcome()))
    testthat::expect_equal(length(clin.res$EW_beta),ncol(mech.res$X),
      label = "Check that the number of betas computed is correct")
})