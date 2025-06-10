test_that("lassoMLE returns expected structure", {
  set.seed(123)
  n <- 50
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)
  lfit <- glmnet::cv.glmnet(X, y, standardize = FALSE, intercept = FALSE)
  expect_error(
    lassoMLE(y, X, lassoFit = lfit,
             optimSteps = 50, sampSteps = 50,
             delay = 1, verbose = FALSE),
    "non-conformable"
  )
})

test_that("invalid method argument fails", {
  set.seed(1)
  n <- 30
  p <- 3
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)
  fit <- glmnet::cv.glmnet(X, y, standardize = FALSE, intercept = FALSE)
  expect_error(
    lassoMLE(y, X, lassoFit = fit, method = "bogus",
             optimSteps = 5, sampSteps = 5,
             delay = 1, verbose = FALSE),
    "Method must be either exact or selected!"
  )
})
