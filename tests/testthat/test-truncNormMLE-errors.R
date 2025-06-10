context("truncNormMLE error handling")

test_that("dimension checks trigger errors", {
  y <- rnorm(3)
  sigma <- matrix(1, 2, 2)
  expect_error(truncNormMLE(y, sigma, threshold = 1, verbose = FALSE),
               "length of y must equal the dimension of sigma!")

  y2 <- rnorm(4)
  sigma2 <- matrix(1, 4, 3)
  expect_error(truncNormMLE(y2, sigma2, threshold = 1, verbose = FALSE),
               "sigma must be a symmetric matrix!")
})

test_that("threshold shape validation", {
  y <- rnorm(2)
  sigma <- diag(2)
  thr <- matrix(1, 3, 2)
  expect_error(truncNormMLE(y, sigma, threshold = thr, verbose = FALSE),
               "threshold must be either a scalar")
})
