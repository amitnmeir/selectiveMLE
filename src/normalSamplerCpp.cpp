#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::plugins(cpp11)]]
#include <progress.hpp>
#include "mleHeader.h"
using namespace Rcpp;
using namespace arma;

/*
 * normalSamplerCpp.cpp
 *
 * Implements a simple Gibbs sampler for multivariate normal vectors
 * under coordinate-wise truncation.  This sampler is used by the
 * `truncNormMLE` routine to approximate the selective MLE when the
 * selection event is described by thresholding.
 */
/**
 * Sample from a normal distribution truncated to the interval
 * `[lower, upper]`.
 */
double sampleBoundedTruncNorm(double mu, double sd, double lower, double upper) {
  double u = runif(1)[0] ;
  double phiB = R::pnorm5(upper, mu, sd, 1, 0) ;
  double phiA = R::pnorm5(lower, mu, sd, 1, 0) ;
  double quantile = u * phiB + phiA * (1 - u) ;
  double sample = R::qnorm(quantile, mu, sd, 1, 0) ;
  return sample ;
}

// [[Rcpp::export]]
/**
 * Gibbs sampler for a multivariate normal vector subject to
 * coordinate-wise truncation.  The `selected` argument indicates which
 * coordinates are truncated on one side only and which are bounded on
 * both sides.  Samples are returned after an initial burn-in period and
 * optional thinning.
 */
NumericVector mvtSampler(NumericVector y,
                         NumericVector mu,
                         IntegerVector selected,
                         NumericMatrix threshold,
                         NumericMatrix precision,
                         int nsamp, int burnin, int trim,
                         bool verbose) {
  int totalIter = burnin + (nsamp - 1) * trim ;
  Progress pb(totalIter, verbose);
  int pdim = y.length() ;
  NumericMatrix samples(nsamp, pdim) ;
  NumericVector samp = clone(y) ;
  double condmean, condsd, pprob, nprob, u ;
  int row = 0;

  for(int i = 0 ; i < (burnin + trim * nsamp) ; i ++) {
    for(int j = 0 ; j < pdim ; j ++) {
      condmean = computeConditionalMean(mu, samp, precision, 1, j) ;
      condsd = 1 / std::sqrt(precision(j, j)) ;
      if(selected[j] == 1) {
        nprob = R::pnorm(threshold(j, 0), condmean, condsd, 1, 1) ;
        pprob = R::pnorm(threshold(j, 1), condmean, condsd, 0, 1) ;
        nprob = 1 / (1 + std::exp(pprob - nprob)) ;
        u = runif(1)[0] ;
        if(u < nprob) {
          samp[j] = sampleUnivTruncNorm(condmean, condsd, threshold(j, 0)) ;
        } else {
          samp[j] = sampleUnivTruncNorm(-condmean, condsd, -threshold(j, 1)) ;
          samp[j] = -samp[j] ;
        }
      } else {
        samp[j] = sampleBoundedTruncNorm(condmean, condsd, threshold(j, 0), threshold(j, 1)) ;
      }
    }

    if(verbose) pb.increment();

    if(i >= (burnin - 1) & (i - burnin + 1) % trim == 0) {
      for(int j = 0 ; j < samples.ncol() ; j++) {
        samples(row, j) = samp[j] ;
      }
      if(++row == nsamp) {
        break ;
      }
    }
  }

  return samples ;
}


