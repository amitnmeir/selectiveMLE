test_that("truncNormMLE basic functionality", {
  y <- c(1.5, rnorm(4))
  sigma <- diag(5)
  res <- truncNormMLE(y, sigma, threshold = 1, maxiter = 20, verbose = FALSE)
  expect_s3_class(res, "truncNormMLE")
  k <- sum(abs(y) > 1)
  expect_length(res$mle, k)
  expect_equal(dim(res$CI), c(k, 2))
})
