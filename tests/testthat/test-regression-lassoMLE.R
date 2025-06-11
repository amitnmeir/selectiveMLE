test_that("lassoMLE replicates with fixed seed", {
  library(glmnet)
  set.seed(123)
  X <- matrix(rnorm(30 * 4), 30, 4)
  beta <- c(1, -0.5, 0, 0)
  y <- as.numeric(X %*% beta + rnorm(30))
  set.seed(123)
  fit1 <- lassoMLE(y, X, optimSteps = 15, sampSteps = 15,
                   delay = 1, verbose = FALSE,
                   multiThread = FALSE)
  set.seed(123)
  fit2 <- lassoMLE(y, X, optimSteps = 15, sampSteps = 15,
                   delay = 1, stepRate = 0.5, verbose = FALSE,
                   multiThread = FALSE)

  expect_equal(fit1$conditionalBeta,
               c(1.11105249, -0.62402436, -0.20422738, 0.06220997),
               tolerance = 1e-6)
  expect_equal(fit2$conditionalBeta,
               c(1.11048631, -0.62360140, -0.20331446, 0.05786993),
               tolerance = 1e-6)

  expect_equal(fit1$wald_CI,
               matrix(c(0.4646590, -0.5930577, -0.4052570, -0.2225908,
                        0.8613642, -0.0344764, 0.1484419, 0.3601090),
                      ncol = 2), tolerance = 1e-6)
  expect_equal(fit2$wald_CI,
               matrix(c(0.4644322, -0.5929839, -0.4049130, -0.2254380,
                        0.8610284, -0.0342704, 0.1491942, 0.3594484),
                      ncol = 2), tolerance = 1e-6)

  expect_equal(fit1$activeSet, fit2$activeSet)
})
