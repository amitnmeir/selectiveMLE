/*
 * lassoSamplerCpp.cpp
 * --------------------
 * This file houses the core C++ routines that power the MCMC based
 * selective maximum likelihood estimator used by the `selectiveMLE`
 * R package.  The implementation relies heavily on Rcpp in order to
 * expose efficient sampling and optimisation code to the R environment.
 *
 * The main entry point is `lassoSampler` which performs a Gibbs style
 * sampler combined with a stochastic optimisation step for estimating
 * the conditional distribution of the lasso regression coefficients.
 * To support this sampler a number of utilities are included:
 *   - truncated normal sampling routines for both extreme and
 *     near-boundary cases;
 *   - convenience methods for shuffling, copying and simple linear
 *     algebra calculations;
 *   - routines for checking and enforcing polyhedral constraints
 *     that arise in the lasso selection event;
 *   - a Metropolis–Hastings step for updating the active signs of the
 *     regression parameters; and
 *   - helpers for computing gradients and bounding the optimisation
 *     trajectory.
 *
 * All functions are written in plain C++ and make extensive use of the
 * Rcpp API for vectorised operations and access to R’s random number
 * generators.  None of the functions are exported directly to R except
 * for `lassoSampler` which is registered via Rcpp attributes.
 */

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]
#include <progress.hpp>
#include <omp.h>
#include <cmath>
#include <algorithm>
using namespace Rcpp;
using namespace arma;

/*
 * lassoSamplerCpp.cpp
 *
 * This file implements the Metropolis--Hastings and Gibbs sampling
 * routines used by the `lassoMLE` R function. The sampler generates
 * draws from the conditional distribution of the refitted regression
 * coefficients given that a particular lasso model was selected.  The
 * samples are then used within a stochastic gradient procedure to
 * approximate the selective maximum likelihood estimator described in
 * Lee et al. (2016).  The code relies heavily on truncated normal
 * proposals and coordinate-wise updates.
 *
 * All functions are designed to be called from R via Rcpp and assume
 * column–major numeric vectors and matrices.
 */

# define M_PI           3.14159265358979323846  /* pi */
const double log2pi = log(2.0 * 3.1415926535897932384);

/**
 * Utility function used primarily for debugging.
 * Prints the elements of a NumericVector to R's
 * console separated by spaces.
 */
void printVec(NumericVector x) {
  for(int i = 0; i < x.length() ; i ++) {
    Rcpp::Rcout<<x[i]<<" ";
  }
  Rcpp::Rcout<<"\n" ;
}

/**
 * Wrapper for R's uniform RNG so that it can be
 * supplied to std::random_shuffle.
 */
inline int randWrapper(const int n) { return floor(unif_rand()*n); }

/**
 * Randomly shuffles the elements of an integer vector
 * in place using R's RNG state.
 */
void randomShuffle(IntegerVector a) {
  std::random_shuffle(a.begin(), a.end(), randWrapper);
}

/**
 * Draws from a truncated normal distribution when the truncation
 * point is far into the tail.  The method uses an exponential
 * rejection sampler which is efficient for extreme values.
 *
 * @param mu        Mean of the normal distribution.
 * @param sd        Standard deviation of the normal distribution.
 * @param threshold Lower bound on the support of the distribution.
 */
double sampleExtreme(double mu, double sd, double threshold) {
  double sign = 1 ;
  double proposal ;
  double alpha ;
  double phi ;

  sign = -1 ;
  mu *= sign ;
  threshold = threshold * sign ;

  // rescaling
  threshold = (threshold - mu) / sd ;
  alpha = (threshold + sqrt(std::pow(threshold, 2) + 4)) / 2 ;

  bool reject = true ;
  int iter = 0;
  while(reject & (iter++ < 10000)) {
    proposal = threshold + R::rexp(1 / alpha) ;
    phi = exp(-std::pow(proposal - alpha, 2) / 2) ;
    if(runif(1)[0] < phi) {
      reject = false ;
    }
  }

  proposal = proposal * sd + mu ;
  return proposal * sign;
}

/**
 * Sample from a univariate truncated normal distribution.  For moderate
 * truncation the sample is obtained by inverting the CDF; for more
 * extreme cases `sampleExtreme` is used as a fall back.
 *
 * @param mu        Mean of the normal distribution.
 * @param sd        Standard deviation of the normal distribution.
 * @param threshold Lower truncation point.
 */
double sampleUnivTruncNorm(double mu, double sd, double threshold) {
  double u = runif(1)[0] ;
  double phiThreshold, sample ;

  // if(isnan(mu)) {
  //   Rcpp::Rcout<<"mu is nan \n" ;
  //   return 0 ;
  // }

  if((std::abs(mu - threshold) / sd) > 3) {
    return sampleExtreme(mu, sd, threshold) ;
  }

  phiThreshold = R::pnorm5(threshold, mu, sd, 1, 0) ;
  sample = R::qnorm5(u * phiThreshold, mu, sd, 1, 0) ;

  int tries = 0 ;
  while(std::isnan(sample) && tries++ < 10) {
    sample = sampleExtreme(mu, sd, threshold) ;
  }

  return sample ;
}

/**
 * Computes the conditional mean of one coordinate of a
 * multivariate normal random vector given the remaining
 * coordinates.  The covariance structure is provided via
 * the matrix `XmX` and the residual variance `yvar`.
 *
 * @param mu    Vector of unconditional means.
 * @param samp  Current state of the sampler.
 * @param XmX   Scaled design cross product matrix.
 * @param yvar  Residual variance (sigma^2).
 * @param index Coordinate for which the mean is computed.
 */
double computeConditionalMean(NumericVector mu,
                              NumericVector samp,
                              const NumericMatrix XmX,
                              double yvar,
                              int index) {
  double result = 0 ;

  for(int j = 0; j < mu.length() ; j ++) {
    if(j != index) {
      result += XmX(index, j) * (samp[j] - mu[j]) / yvar ;
    }
  }

  result = result / XmX(index, index) * yvar ;
  result = mu[index] + result ;
  return result ;
}

/**
 * Check whether every element of `samp` lies between the
 * corresponding lower `l` and upper `u` bounds.  Returns 1 if the
 * condition holds and 0 otherwise.  Used to verify whether the current
 * sample satisfies the zero-coordinate constraints.
 */
int innerZeroCondition(NumericVector samp,
                       const NumericVector l,
                       const NumericVector u) {
  int count = 1;
  double x ;

  for(int j = 0; j < samp.length() ; j++) {
    x = samp[j] ;
    if(x > u[j] | x < l[j]) {
      count-- ;
      break ;
    }
  }

  return count;
}

/**
 * Determine whether the current sample of coefficients satisfies the
 * one-step sampling constraints encoded in `threshold`.  Returns 1 if
 * all coordinates are within the allowed region and 0 otherwise.
 */
int innerCheckOneCondition(NumericVector samp,
                           NumericVector threshold) {
  double x ;
  for(int i = 0; i < samp.length() ; i++) {
    x = - std::abs(samp[i]) ;
    if(x > threshold[i]) {
      return 0 ;
    }
  }

  return 1;
}

/**
 * Compute the lower and upper bounds defining the selection region for
 * the coordinates that were not selected by the lasso.  The matrix
 * `u0mat` corresponds to the design matrix of inactive variables and
 * `signs` gives the current sign configuration of the active
 * coefficients.
 */
void computeZeroThresholds(const NumericMatrix u0mat,
                           NumericVector signs,
                           NumericVector l,
                           NumericVector u) {
  if(u0mat.ncol() != static_cast<int>(signs.size())) {
    #pragma omp parallel for default(none) shared(u, l, u0mat)
    for(int i = 0; i < u0mat.nrow(); ++i) {
      u[i] = 1.0;
      l[i] = -1.0;
    }
    return;
  }
  arma::mat A(u0mat.begin(), u0mat.nrow(), u0mat.ncol());
  arma::colvec s(signs.begin(), signs.size(), false);
  arma::colvec prod = A * s;
  #pragma omp parallel for default(none) shared(u, l, prod, A)
  for(int i = 0; i < A.n_rows; ++i) {
    u[i] = 1.0 - prod[i];
    l[i] = -1.0 - prod[i];
  }
}

/**
 * Compute the truncation thresholds used when sampling the active
 * coefficients.  The thresholds depend on the current sign vector and
 * the inverse design matrix `XmXinv` scaled by the lasso penalty
 * `lambda`.
 */
void computeOneThreshold(NumericVector signs,
                         double lambda,
                         const NumericMatrix &XmXinv,
                         NumericVector u) {
  arma::mat inv(XmXinv.begin(), XmXinv.nrow(), XmXinv.ncol());
  arma::colvec s(signs.begin(), signs.size(), false);
  arma::colvec res = lambda * inv * s;
  std::copy(res.begin(), res.end(), u.begin());
}

/**
 * Compute the change in the truncation threshold that would result from
 * flipping the sign of a single coordinate.  This is used when
 * proposing sign changes in the Metropolis--Hastings step.
 */
double computeDiffThreshold(NumericVector signs,
                            double lambda,
                            const NumericMatrix &XmXinv,
                            int coordinate) {
  arma::mat inv(XmXinv.begin(), XmXinv.nrow(), XmXinv.ncol());
  arma::colvec s(signs.begin(), signs.size(), false);
  double rowDot = arma::dot(inv.row(coordinate), s);
  double result = rowDot - 2.0 * inv(coordinate, coordinate) * signs[coordinate];
  return result * lambda;
}

/** Copy the contents of one numeric vector to another. */
void copyVector(NumericVector &to, const NumericVector &from) {
  std::copy(from.begin(), from.end(), to.begin());
}

/**
 * Gibbs updates for the active coefficients.
 * Each coordinate is sampled from its truncated
 * normal full conditional distribution.
 *
 * @param samp      Current coefficient vector (updated in place).
 * @param signs     Sign vector of the active coefficients.
 * @param u         Thresholds derived from the lasso penalty.
 * @param mean      Conditional mean of the regression coefficients.
 * @param XmX       Cross product matrix of the design.
 * @param condSigma Conditional variances for each coefficient.
 * @param ysigsq    Residual variance.
 * @param burnin    Unused but kept for symmetry with R version.
 * @param maxiter   Number of Gibbs sweeps to perform.
 */
void aOneSampler(NumericVector &samp,
                 const NumericVector &signs,
                 const NumericVector &u,
                 const NumericVector &mean,
                 const NumericMatrix &XmX,
                 const NumericVector &condSigma,
                 double ysigsq,
                 int burnin, int maxiter) {
  double condMean ;
  double condSD ;
  int i, j ;

  int iterations = 0 ;
  while(iterations++ < maxiter) {
    for(int j = 0; j < samp.length() ; j++) {
      condMean = computeConditionalMean(mean, samp, XmX, ysigsq, j) ;
      condSD = std::sqrt(condSigma[j]) ;
      samp[j] = sampleUnivTruncNorm(-signs[j] * condMean, condSD, -signs[j] * u[j]) ;
      samp[j] *= -signs[j] ;
    }
  }
}

/**
 * Block Gibbs sampler for the components of Ay that correspond
 * to inactive variables (the so called "zero" part).  Sampling
 * is performed under box constraints using a truncated Gaussian
 * proposal obtained via sequential univariate updates.
 */
void zeroSampler(NumericVector &zsamp,
                 NumericVector &zeroSamp,
                 NumericVector zeroMean,
                 NumericMatrix sqrtMat,
                 NumericVector lzero,
                 NumericVector uzero,
                 int miniters, int maxiters) {
  double eps = 0.00000000001 ;

  double currentLower, currentUpper ;
  double threshold ;
  double adjustment ;
  double x, sigmaij ;
  int i, j ;

  for(int i = 0; i < zsamp.length() ; i++) zsamp[i] = 0 ;
  for(int i = 0; i < zeroMean.length() ; i++) zeroSamp[i] = zeroMean[i] ;

  for(int k = 0 ; k < maxiters ; k++) {
    for(int i = 0 ; i < zsamp.length() ; i ++) {
      // Finding constraints for next sample
      currentLower = -10000000000 ;
      currentUpper = 10000000000 ;
      for(int j = 0; j < zeroSamp.length() ; j ++) {
        sigmaij = sqrtMat(j, i) ;
        if(std::abs(sigmaij) < eps) continue ;
        adjustment = - zeroSamp[j] + sigmaij * zsamp[i] ;
        if(sigmaij > 0) {
          currentLower = std::max(currentLower, (lzero[j] + adjustment) / sigmaij) ;
          currentUpper = std::min(currentUpper, (uzero[j] + adjustment) / sigmaij);
        } else {
          currentLower = std::max(currentLower, (uzero[j] + adjustment) / sigmaij) ;
          currentUpper = std::min(currentUpper, (lzero[j] + adjustment) / sigmaij);
        }
      }

      // truncated sampling
      x = R::pnorm5(currentUpper, 0, 1, 1, 0) - R::pnorm5(currentLower, 0, 1, 1, 0);
      x *= runif(1)[0] ;
      x += R::pnorm5(currentLower, 0, 1, 1, 0) ;
      x = R::qnorm5(x, 0, 1, 1, 0) ;
      for(int j = 0; j < sqrtMat.nrow() ; j++) {
        sigmaij = sqrtMat(j, i) ;
        zeroSamp[j] += x * sigmaij - sigmaij * zsamp[i];
      }
      zsamp[i] = x ;
    }
    if(k > miniters) {
      if(innerZeroCondition(zeroSamp, lzero, uzero) == 1) {
        return ;
      }
    }
  }
}

/**
 * Utility to obtain the sign of a scalar with a small epsilon
 * region around zero treated as exactly zero.
 */
double sign(double x) {
  double eps = 0.00000001 ;
  if(std::abs(x) < eps) {
    return 0.0 ;
  } else if(x < 0) {
    return -1.0 ;
  } else {
    return 1.0 ;
  }
}

/**
 * Compute a stochastic gradient step for the selective likelihood.  The
 * gradient is scaled by a decreasing step size controlled by
 * `stepCoef`, `stepRate` and `delay` and is clipped to
 * `gradientBound` to ensure stability.
 */
void computeGradient(NumericVector gradient, NumericVector Xy,
                     NumericVector samp, NumericMatrix XmX,
                     double stepCoef, double stepRate,
                     double gradientBound,
                     int iter, int delay) {
  stepCoef = stepCoef / std::pow(std::max(1, iter + 1 - delay), stepRate) ;
#pragma omp parallel for default(none) shared(gradient, Xy, samp, XmX, stepCoef, gradientBound)
  for(int i = 0 ; i < gradient.length() ; i++) {
    double g = Xy[i];
    for(int j = 0 ; j < XmX.ncol() ; j++) {
      g -= XmX(i, j) * samp[j] ;
    }
    g *= stepCoef ;
    if(std::abs(g) > gradientBound) {
      g = gradientBound * sign(g) ;
    }
    gradient[i] = g;
  }
}

/**
 * Compute a running mean of the estimated coefficients over the last
 * 300 iterations of the optimization.  The result is stored in
 * `betaOut`.
 */
void computeConditionalBeta(NumericVector &betaOut,
                            const NumericMatrix &estimateMat,
                            int iter) {
  int meanStart = std::max(0, iter - 300);
  #pragma omp parallel for default(none) shared(betaOut, estimateMat, iter, meanStart)
  for(int i = 0; i < estimateMat.ncol(); i++) {
    double sum = 0.0;
    int denominator = 0;
    for(int j = meanStart; j <= iter; j++) {
      denominator++;
      sum += estimateMat(j, i);
    }
    betaOut[i] = sum / denominator;
  }
}

/**
 * Ensures that the optimiser stays within the box implied by the
 * naive lasso solution.  Coefficients are clipped so that they do
 * not cross zero or exceed the naive estimates in magnitude.
 */
void boundBeta(NumericVector &estimate, NumericVector naive) {
  #pragma omp parallel for default(none) shared(estimate, naive)
  for(int i = 0 ; i < estimate.length() ; i++) {
    if(naive[i] < 0) {
      if(estimate[i] < naive[i]) {
        estimate[i] = naive[i] ;
      } else if(estimate[i] > 0) {
        estimate[i] = 0;
      }
    } else {
      if(estimate[i] > naive[i]) {
        estimate[i] = naive[i] ;
      } else if(estimate[i] < 0) {
        estimate[i] = 0 ;
      }
    }
  }
}

/**
 * Evaluate the log-density of a multivariate normal distribution
 * with covariance proportional to the inverse of `XmX`.
 * This is used in the Metropolis-Hastings ratio for sign updates.
 */
double mvtLogDens(NumericVector samp, NumericVector mean,
                  NumericMatrix XmX, double ysigsq) {
  NumericVector diff(samp.length());
  #pragma omp parallel for default(none) shared(diff, samp, mean)
  for(int i = 0 ; i < diff.length() ; i++) {
    diff[i] = samp[i] - mean[i];
  }

  double result = 0.0;
  #pragma omp parallel for reduction(+:result) default(none) shared(diff, XmX)
  for(int i = 0 ; i < diff.length() ; i++) {
    double inner = 0.0;
    for(int j = 0 ; j < diff.length() ; j++) {
      inner += XmX(i, j) * diff[j];
    }
    result += inner * diff[i];
  }

  return -0.5 / ysigsq * result;
}

/**
 * Metropolis--Hastings update that proposes sign changes for the active
 * coefficients and possibly resamples all coefficients if necessary.
 * This step is responsible for exploring the truncated distribution
 * defined by the lasso selection event.
 */
void mhSampler(NumericVector samp, NumericVector oldSamp,
               NumericVector signs, NumericVector newsigns,
               NumericVector mean, NumericMatrix &sigma,
               NumericVector condSigma,
               double lambda,
               NumericVector u, NumericVector newuone,
               NumericMatrix &XmXinv,
               NumericMatrix &XmX, double ysigsq,
               NumericVector lZero, NumericVector uZero,
               NumericVector newlzero, NumericVector newuzero,
               NumericMatrix &u0mat, NumericVector zeroSamp,
               int maxiter, IntegerVector order, NumericVector initEst,
               bool methodExact) {
  int j ;
  double condSD, condMean ;
  double thres, psame, pdiff ;
  int zeroCondition, oneCondition ;
  double newsamp ;
  double forwardDens, reverseDens, sd;
  double oldMVTdens, newMVTdens, mhRatio ;

  for(int iter = 0 ; iter < maxiter ; iter++) {
    randomShuffle(order) ;
    for(int i = 0 ; i < samp.length() ; i ++) {
      // Deciding whether to change signs
      j = order[i] ;
      if(std::abs(mean[j]) > std::abs(initEst[j]) &
         signs[j] == sign(initEst[j])) continue ;
      newsigns[j] *= - 1 ;
      condMean = computeConditionalMean(mean, samp, XmX, ysigsq, j) ;
      condSD = std::sqrt(condSigma[j]) ;
      thres  = computeDiffThreshold(newsigns, lambda, XmXinv, j) ;
      psame = R::pnorm(-signs[j] * u[j], -signs[j] * condMean, condSD, 1, 1) ;
      pdiff = R::pnorm(signs[j] * thres, signs[j] * condMean, condSD, 1, 1) ;
      pdiff = 1.0 / (1.0 + std::exp(psame - pdiff)) ;
      if(runif(1)[0] > pdiff) {
        newsamp = sampleUnivTruncNorm(-signs[j] * condMean, condSD, - signs[j] * u[j]) ;
        samp[j] = -signs[j] * newsamp ;
        oldSamp[j] = samp[j] ;
        newsigns[j] *= -1 ;
        continue ;
      }

      // check zero condition
      computeZeroThresholds(u0mat, newsigns, newlzero, newuzero) ;
      if(methodExact) {
        zeroCondition = innerZeroCondition(zeroSamp, newlzero, newuzero) ;
      }

      if(methodExact & (zeroCondition == 0)) {
        newsamp = sampleUnivTruncNorm(-signs[j] * condMean, condSD, - signs[j] * u[j]) ;
        samp[j] = -signs[j] * newsamp ;
        oldSamp[j] = samp[j] ;
        newsigns[j] *= -1 ;
        continue ;
      }

      // Check if other coordinates need to be sampled
      computeOneThreshold(newsigns, lambda, XmXinv, newuone) ;
      oneCondition = innerCheckOneCondition(samp, newuone) ;
      if(oneCondition == 1) { // if no, just sample the one coordinate
        signs[j] *= - 1;
        copyVector(u, newuone) ;
        newsamp = sampleUnivTruncNorm(-signs[j] * condMean, condSD, -signs[j] * u[j]) ;
        newsamp *= -signs[j] ;
        samp[j] = newsamp ;
        oldSamp[j] = newsamp ;
        if(methodExact) {
          copyVector(lZero, newlzero) ;
          copyVector(uZero, newuzero) ;
        }
        continue ;
      }

      // If yes, then do MH step
      forwardDens = 0;
      reverseDens = 0;
      for(int k = 0; k < samp.length() ; k++) {
        if(k == j) {
          newsamp = sampleUnivTruncNorm(signs[j] * condMean, condSD, signs[j] * thres) ;
          newsamp *= signs[j] ;
          samp[j] = newsamp ;
          forwardDens += R::dnorm(newsamp, condMean, condSD, 1) ;
        } else {
          sd = std::sqrt(sigma(k, k)) ;
          newsamp = sampleUnivTruncNorm(-signs[k] * samp[k], sd, -signs[k] * newuone[j]) ;
          newsamp *= -signs[k] ;
          samp[k] = newsamp ;
          forwardDens += R::dnorm(newsamp, oldSamp[k], sd, 1) ;
          forwardDens -= R::pnorm(-signs[k] * newuone[k], -signs[k] * oldSamp[k], sd, 1, 1) ;
          reverseDens += R::dnorm(oldSamp[k], samp[k], sd, 1) ;
          reverseDens -= R::pnorm(-signs[k] * u[k], -signs[k] * samp[k], sd, 1, 1) ;
        }
      }
      condMean = computeConditionalMean(mean, samp, XmX, ysigsq, j) ;
      reverseDens += R::dnorm(oldSamp[j], condMean, condSD, 1) ;
      newMVTdens = mvtLogDens(samp, mean, XmX, ysigsq) ;
      oldMVTdens = mvtLogDens(oldSamp, mean, XmX, ysigsq) ;

     // Do we boldly go where no sampler has gone before?
      mhRatio = newMVTdens - oldMVTdens + reverseDens - forwardDens ;
      mhRatio = std::exp(mhRatio) ;
      //Rcpp::Rcout<<" MH "<<mhRatio ;
      if(runif(1)[0] < mhRatio) {
        signs[j] *= -1 ;
        copyVector(oldSamp, samp) ;
        copyVector(u, newuone) ;
        if(methodExact) {
          copyVector(lZero, newlzero) ;
          copyVector(uZero, newuzero) ;
        }
      } else {
        newsigns[j] *= -1 ;
        copyVector(samp, oldSamp) ;
      }
    }
  }
}

/**
 * Main driver routine used by `lassoMLE`.  Runs a block Gibbs sampler
 * combined with stochastic gradient updates in order to approximate the
 * selective MLE for a lasso-selected model.
 */
// [[Rcpp::export]]
/**
 * Main workhorse implementing the selective MLE sampler for the
 * lasso.  The function alternates between sampling the active set
 * coefficients and performing stochastic gradient updates.  It is
 * exposed to R via Rcpp and is not intended to be called directly
 * by users.
 *
 * @param initEst  Initial estimate of the regression coefficients.
 * @param initSamp Starting point for the sampler.
 * @param oneCov   Conditional covariance matrix for the active set.
 * @param XmX      Cross product of the design matrix.
 * @param XmXinv   Inverse of XmX restricted to the active set.
 * @param condSigma Conditional variances for each coefficient.
 * @param lambda   Lasso penalty value.
 * @param ysigsq   Residual variance.
 * @param zeroMean Mean vector for the zero part of Ay.
 * @param sqrtZero Square root of the covariance for the zero part.
 * @param u0mat    Constraint matrix defining the selection event.
 * @param n,p      Dimensions of the design matrix.
 * @param nsamp    Number of Gibbs samples per outer iteration.
 * @param burnin   Length of the burn-in period.
 * @param Xy       Cross product of X and y.
 * @param estimateMat Matrix to store estimates at each iteration.
 * @param sampMat  Matrix to store sampler draws at each iteration.
 * @param delay,stepRate,stepCoef Parameters controlling the optimiser.
 * @param gradientBound Bound on the gradient magnitude.
 * @param assumeConvergence Iteration at which optimisation stops.
 * @param naive    Naive lasso estimate used for bounding.
 * @param methodExact Whether to perform exact polyhedral sampling.
 * @param verbose  If true, progress is printed to the console.
 */
NumericVector lassoSampler(const NumericVector initEst,
                    const NumericVector initSamp,
                    NumericMatrix oneCov, NumericMatrix XmX,
                    NumericMatrix XmXinv,
                    NumericVector condSigma,
                    double lambda, double ysigsq,
                    NumericVector zeroMean, NumericMatrix sqrtZero,
                    NumericMatrix u0mat,
                    int n, int p,
                    int nsamp, int burnin,
                    NumericVector Xy,
                    NumericMatrix &estimateMat, NumericMatrix &sampMat,
                    int delay, double stepRate, double stepCoef,
                    double gradientBound, int assumeConvergence,
                    NumericVector naive, bool methodExact, bool verbose) {
  // Initializing Sampling order
  IntegerVector order = IntegerVector(initEst.length()) ;
  for(int i = 0; i < order.length() ; i++) order[i] = i ;

  // initalizing signs
  NumericVector samp = clone(initSamp) ;
  NumericVector oldsamp = clone(initSamp) ;
  NumericVector estimate = clone(initEst) ;
  NumericVector signs = NumericVector(samp.length()) ;
  for(int i = 0; i < signs.length() ; i++) {
    signs[i] = -1 ;
    if(samp[i] > 0) signs[i] = 1 ;
  }
  NumericVector newsigns = clone(signs) ;

  // Initializing zero sample
  int k = samp.length() ;
  NumericVector zsamp, zeroSamp, lzero, uzero, newlzero, newuzero ;
  if(methodExact) {
    int zeroRank = std::min(p - k, n - k) ;
    zsamp = NumericVector(zeroRank) ;
    zeroSamp = NumericVector(p - k) ;
    lzero = NumericVector(p - k) ;
    uzero = NumericVector(p - k) ;
    newlzero = NumericVector(p - k) ;
    newuzero = NumericVector(p - k) ;
    computeZeroThresholds(u0mat, signs, lzero, uzero) ;
  }

  // initializing one sampler
  double currentThreshold ;
  NumericVector uone = NumericVector(samp.length()) ;
  NumericVector newuone = NumericVector(samp.length()) ;
  computeOneThreshold(signs, lambda, XmXinv, uone) ;

  // initializing optimizer
  NumericVector gradient = NumericVector(samp.length()) ;
  Progress pb(sampMat.nrow(), verbose);
  int frac = std::max(1, sampMat.nrow() / 20);

  for(int optimIter = 0 ; optimIter < sampMat.nrow() ; optimIter ++) {
    // SAMPLING STARTS
    for(int sampIter = 0 ; sampIter < nsamp ; sampIter++){
      // sampling A0y
      if(methodExact) {
        zeroSampler(zsamp, zeroSamp,
                    zeroMean, sqrtZero,
                    lzero, uzero, 5, 40) ;
      }

      // Sampling Signs
      mhSampler(samp, oldsamp, signs, newsigns,
                estimate, oneCov, condSigma,
                lambda, uone, newuone,
                XmXinv, XmX, ysigsq,
                lzero, uzero, newlzero, newuzero,
                u0mat, zeroSamp,
                1, order, initEst, methodExact) ;

      // Sampling regression Coefs
      aOneSampler(samp, signs, uone, estimate, XmX, condSigma, ysigsq, 10, 10) ;
      copyVector(oldsamp, samp) ;
    }

    // Rcpp::Rcout<<optimIter<<" ";
    // printVec(samp) ;

    // Checking that Sample is reasonable
    for(int l = 0 ; l < samp.length() ; l++) {
      if(std::abs((naive[l] - samp[l]) / std::sqrt(oneCov(l, l))) > 10) {
        for(int i = 0 ; i < samp.length() ; i ++) {
          samp[i] = naive[i] ;
          signs[i] = sign(naive[i]) ;
        }
        copyVector(oldsamp, samp) ;
        copyVector(newsigns, signs) ;
        computeOneThreshold(signs, lambda, XmXinv, uone) ;
        if(methodExact) {
          computeZeroThresholds(u0mat, signs, lzero, uzero) ;
        }
        break ;
      }
    }

    // SAMPLING ENDS

    // OPTIMIZATION STARTS
    // computing gradient and updating estimate
    if(optimIter < assumeConvergence) {
      // Beta Update
      computeGradient(gradient, Xy, samp, XmX,
                      stepCoef, stepRate, gradientBound,
                      optimIter, delay) ;
      for(int i = 0 ; i < gradient.length() ; i++) {
        estimate[i] += gradient[i] ;
      }
      boundBeta(estimate, naive) ;
    }

    // `reporting' sample and estimate
    for(int i = 0 ; i < sampMat.ncol() ; i++) {
      sampMat(optimIter, i) = samp[i] ;
      estimateMat(optimIter, i) = estimate[i] ;
    }

    if(optimIter == (assumeConvergence - 1)) {
      computeConditionalBeta(estimate, estimateMat, optimIter) ;
      burnin *= 2 ;
    }

    // OPTIMIZATION ENDS
    if (verbose) pb.increment();
    if (verbose && (optimIter + 1) % frac == 0) {
      for(int i = 0 ; i < samp.length() ; i ++) {
        samp[i] = naive[i] ;
        signs[i] = sign(naive[i]) ;
      }
      copyVector(oldsamp, samp) ;
      copyVector(newsigns, signs) ;
      computeOneThreshold(signs, lambda, XmXinv, uone) ;
      if(methodExact) {
        computeZeroThresholds(u0mat, signs, lzero, uzero) ;
      }
    }
  }

  return samp;
}
