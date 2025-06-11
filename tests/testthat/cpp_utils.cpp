#include <Rcpp.h>
// [[Rcpp::plugins(openmp)]]
#include <cmath>
#include <algorithm>
#include <omp.h>
using namespace Rcpp;

/*
 * Wrappers for sampler utility functions used in the C++ implementation of
 * the lasso sampler. The original functions are defined in
 * `src/lassoSamplerCpp.cpp` but are not exported to R. For the unit tests we
 * replicate their bodies here and expose thin wrappers with `Rcpp::export` so
 * the helpers can be called directly from R. This keeps the package interface
 * unchanged while still allowing the low level algorithms to be verified.
 */

double computeConditionalMean_cpp(NumericVector mu,
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

int innerZeroCondition_cpp(NumericVector samp,
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

int innerCheckOneCondition_cpp(NumericVector samp,
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

void computeZeroThresholds_cpp(const NumericMatrix u0mat,
                               NumericVector signs,
                               NumericVector l,
                               NumericVector u) {
  #pragma omp parallel for default(none) shared(u, l, signs, u0mat)
  for(int i = 0; i < u0mat.nrow() ; i ++) {
    u[i] = 1;
    l[i] = -1;
    for(int j = 0; j < signs.length() ; j ++) {
      u[i] -= u0mat(i, j) * signs[j] ;
      l[i] -= u0mat(i, j) * signs[j] ;
    }
  }
}

void computeOneThreshold_cpp(NumericVector signs,
                             double lambda,
                             const NumericMatrix XmXinv,
                             NumericVector u) {
  #pragma omp parallel for default(none) shared(u, signs, XmXinv, lambda)
  for(int i = 0 ; i < XmXinv.nrow() ; i++) {
    double val = 0 ;
    for(int j = 0 ; j < XmXinv.ncol() ; j++) {
      val += XmXinv(i, j) * signs[j] ;
    }
    u[i] = val * lambda ;
  }
}

double computeDiffThreshold_cpp(NumericVector signs,
                                double lambda,
                                const NumericMatrix XmXinv,
                                int coordinate) {
  double result = 0;
  for(int i = 0 ; i < XmXinv.ncol() ; i++) {
    if(i == coordinate) {
      result -= XmXinv(coordinate, i) * signs[i] ;
    } else {
      result += XmXinv(coordinate, i) * signs[i] ;
    }
  }
  return result * lambda ;
}

double sign_cpp(double x) {
  double eps = 0.00000001 ;
  if(std::abs(x) < eps) {
    return 0.0 ;
  } else if(x < 0) {
    return -1.0 ;
  } else {
    return 1.0 ;
  }
}

void boundBeta_cpp(NumericVector estimate, NumericVector naive) {
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

void computeGradient_cpp(NumericVector gradient, NumericVector Xy,
                         NumericVector samp, NumericMatrix XmX,
                         double stepCoef, double stepRate,
                         double gradientBound,
                         int iter, int delay) {
  stepCoef = stepCoef / std::pow(std::max(1, iter + 1 - delay), stepRate);
#pragma omp parallel for default(none) shared(gradient, Xy, samp, XmX, stepCoef, gradientBound)
  for(int i = 0 ; i < gradient.length() ; i++) {
    double g = Xy[i];
    for(int j = 0 ; j < XmX.ncol() ; j++) {
      g -= XmX(i, j) * samp[j];
    }
    g *= stepCoef;
    if(std::abs(g) > gradientBound) {
      g = gradientBound * sign_cpp(g);
    }
    gradient[i] = g;
  }
}

// Exported wrappers ---------------------------------------------------------

/* Compute conditional mean (test wrapper)
 *
 * This function simply forwards to `computeConditionalMean_cpp` which
 * mirrors the internal sampler helper.  Exposed for unit tests only.
 */
// [[Rcpp::export]]
double test_computeConditionalMean(NumericVector mu, NumericVector samp,
                                   NumericMatrix XmX, double yvar, int index) {
  return computeConditionalMean_cpp(mu, samp, XmX, yvar, index);
}

/* Check zero-condition bounds (test wrapper)
 * Exposes `innerZeroCondition_cpp` to R so tests can verify the
 * selection-region logic.
 */
// [[Rcpp::export]]
int test_innerZeroCondition(NumericVector samp, NumericVector l,
                            NumericVector u) {
  return innerZeroCondition_cpp(samp, l, u);
}

/* Check one-condition thresholds (test wrapper) */
// [[Rcpp::export]]
int test_innerCheckOneCondition(NumericVector samp, NumericVector thresh) {
  return innerCheckOneCondition_cpp(samp, thresh);
}

/* Compute zero-coordinate thresholds (test wrapper) */
// [[Rcpp::export]]
NumericMatrix test_computeZeroThresholds(NumericMatrix u0mat,
                                         NumericVector signs) {
  NumericVector l(u0mat.nrow()), u(u0mat.nrow());
  computeZeroThresholds_cpp(u0mat, signs, l, u);
  NumericMatrix out(u0mat.nrow(),2);
  for(int i=0;i<u0mat.nrow();i++){
    out(i,0)=l[i];
    out(i,1)=u[i];
  }
  return out;
}

/* Compute one-coordinate thresholds (test wrapper) */
// [[Rcpp::export]]
NumericVector test_computeOneThreshold(NumericVector signs, double lambda,
                                       NumericMatrix XmXinv) {
  NumericVector u(signs.length());
  computeOneThreshold_cpp(signs, lambda, XmXinv, u);
  return u;
}

/* Compute threshold difference for a coordinate (test wrapper) */
// [[Rcpp::export]]
double test_computeDiffThreshold(NumericVector signs, double lambda,
                                 NumericMatrix XmXinv, int coordinate) {
  return computeDiffThreshold_cpp(signs, lambda, XmXinv, coordinate);
}

/* Sign helper used in gradient updates (test wrapper) */
// [[Rcpp::export]]
double test_sign(double x) {
  return sign_cpp(x);
}

/* Bound estimated coefficients by naive solution (test wrapper) */
// [[Rcpp::export]]
NumericVector test_boundBeta(NumericVector estimate, NumericVector naive) {
  NumericVector est = clone(estimate);
  boundBeta_cpp(est, naive);
  return est;
}

/* Compute gradient step for sampler optimisation (test wrapper) */
// [[Rcpp::export]]
NumericVector test_computeGradient(NumericVector gradient, NumericVector Xy,
                                   NumericVector samp, NumericMatrix XmX,
                                   double stepCoef, double stepRate,
                                   double gradientBound, int iter, int delay) {
  NumericVector grad = clone(gradient);
  computeGradient_cpp(grad, Xy, samp, XmX, stepCoef, stepRate,
                      gradientBound, iter, delay);
  return grad;
}

