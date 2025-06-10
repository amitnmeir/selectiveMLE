# Tests for the sampler utility functions implemented in C++.
# The utilities themselves are duplicated in `cpp_utils.cpp` and compiled
# via `helper-cpp-utils.R`.  Exporting them this way keeps the package API
# untouched while allowing `testthat` to verify the low level calculations.
context('C++ sampler utilities')

# compile helpers via helper-cpp-utils.R (already sourced)

test_that('computeConditionalMean works', {
  mu <- c(0, 0)
  samp <- c(1, 2)
  XmX <- diag(2)
  expect_equal(test_computeConditionalMean(mu, samp, XmX, 1, 1), 0)
})

test_that('innerZeroCondition identifies bounds', {
  samp <- c(0.5, -0.2)
  l <- c(-1, -1)
  u <- c(1, 1)
  expect_equal(test_innerZeroCondition(samp, l, u), 1)
  samp2 <- c(1.5, 0)
  expect_equal(test_innerZeroCondition(samp2, l, u), 0)
})

test_that('innerCheckOneCondition behaves correctly', {
  samp <- c(0.2, -0.2)
  thr <- c(-0.1, -0.1)
  expect_equal(test_innerCheckOneCondition(samp, thr), 1)
  thr2 <- c(-0.3, 0.1)
  expect_equal(test_innerCheckOneCondition(samp, thr2), 0)
})

test_that('threshold computations are correct', {
  u0mat <- diag(2)
  signs <- c(1, -1)
  res <- test_computeZeroThresholds(u0mat, signs)
  expect_equal(res[,1], c(-2, 0))
  expect_equal(res[,2], c(0, 2))
  XmXinv <- diag(2)
  expect_equal(test_computeOneThreshold(signs, 2, XmXinv), c(2, -2))
  expect_equal(test_computeDiffThreshold(signs, 2, XmXinv, 1), 2)
})

test_that('sign and bounding utilities work', {
  expect_equal(test_sign(0), 0)
  expect_equal(test_sign(1e-9), 0)
  expect_equal(test_sign(-1e-7), -1)
  expect_equal(test_sign(2), 1)

  est <- c(1, -0.5, 0.3)
  naive <- c(0.8, -1, 0.1)
  bounded <- test_boundBeta(est, naive)
  expect_equal(bounded, c(0.8, -0.5, 0.1))
})

test_that('computeGradient produces expected output', {
  grad <- c(0, 0)
  Xy <- c(1, 2)
  samp <- c(0.2, -0.1)
  XmX <- diag(2)
  res <- test_computeGradient(grad, Xy, samp, XmX, 0.1, 0, 100, 0, 0)
  expect_equal(res, c(0.08, 0.21))
})
