# Utility functions for lassoMLE and truncNormMLE

validate_lasso_inputs <- function(y, X, method, lambda) {
  if(!is.numeric(y)) {
    stop("y must be numeric")
  }
  if(!is.matrix(X) || !is.numeric(X)) {
    stop("X must be a numeric matrix")
  }
  if(nrow(X) != length(y)) {
    stop("X and y must have the same number of rows")
  }
  match.arg(method, c("exact", "selected"))
  if(!(is.numeric(lambda) || lambda %in% c("lambda.min", "lambda.1se"))) {
    stop("lambda must be numeric or one of 'lambda.min' or 'lambda.1se'")
  }
  invisible(NULL)
}

standardize_xy <- function(y, X) {
  sdy <- sd(y)
  y_std <- (y / sdy) - mean(y / sdy)
  meanX <- colMeans(X)
  sdX <- apply(X, 2, sd)
  X_std <- sweep(X, 2, meanX, "-")
  X_std <- sweep(X_std, 2, sdX, "/")
  list(y = y_std, X = X_std, sdy = sdy, meanX = meanX, sdX = sdX)
}

compute_lambda <- function(lassoFit, lambda, n) {
  if(lambda == "lambda.min") {
    lassoFit$lambda.min * n
  } else if(lambda == "lambda.1se") {
    lassoFit$lambda.1se * n
  } else {
    as.numeric(lambda)
  }
}

compute_wald <- function(betaSample, XmX, ysig, assumeConvergence, n, conditionalBeta) {
  ysamp <- betaSample %*% XmX / ysig^2
  idx <- seq(assumeConvergence, nrow(betaSample))
  forQuantiles <- ysamp[idx, , drop = FALSE] / sqrt(n)
  center <- colMeans(forQuantiles)
  variance <- var(forQuantiles)
  coefVar <- solve(variance)
  centered <- sweep(forQuantiles, 2, center)
  scaled <- centered %*% coefVar
  q <- t(apply(scaled, 2, function(x) quantile(x, c(0.025, 0.975))))
  wald <- cbind(conditionalBeta - q[,2] / sqrt(n),
                conditionalBeta - q[,1] / sqrt(n))
  list(wald = wald, coefVar = coefVar, variance = variance)
}

validate_trunc_inputs <- function(y, sigma, threshold) {
  if(!is.numeric(y)) {
    stop("y must be numeric")
  }
  if(ncol(sigma) != nrow(sigma)) {
    stop("sigma must be a symmetric matrix!")
  }
  if(length(y) != ncol(sigma)) {
    stop("length of y must equal the dimension of sigma!")
  }
  invisible(NULL)
}

expand_threshold <- function(threshold, p) {
  if(length(threshold) == 1) {
    th <- abs(threshold)
    cbind(rep(-th, p), rep(th, p))
  } else if(length(threshold) == p) {
    cbind(-abs(threshold), abs(threshold))
  } else if(is.matrix(threshold) && all(dim(threshold) == c(p, 2))) {
    threshold
  } else {
    stop("threshold must be either a scalar, a vector of size length(y) or a matrix of size length(y) X 2")
  }
}
