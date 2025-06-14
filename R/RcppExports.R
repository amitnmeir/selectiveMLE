# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

lassoSampler <- function(initEst, initSamp, oneCov, XmX, XmXinv, condSigma, lambda, ysigsq, zeroMean, sqrtZero, u0mat, n, p, nsamp, burnin, Xy, estimateMat, sampMat, delay, stepRate, stepCoef, gradientBound, assumeConvergence, naive, methodExact, verbose) {
    .Call(`_selectiveMLE_lassoSampler`, initEst, initSamp, oneCov, XmX, XmXinv, condSigma, lambda, ysigsq, zeroMean, sqrtZero, u0mat, n, p, nsamp, burnin, Xy, estimateMat, sampMat, delay, stepRate, stepCoef, gradientBound, assumeConvergence, naive, methodExact, verbose)
}

mvtSampler <- function(y, mu, selected, threshold, precision, nsamp, burnin, trim, verbose) {
    .Call(`_selectiveMLE_mvtSampler`, y, mu, selected, threshold, precision, nsamp, burnin, trim, verbose)
}

