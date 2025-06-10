#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

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
  for(int i = 0 ; i < XmXinv.nrow() ; i++) {
    u[i] = 0 ;
    for(int j = 0 ; j < XmXinv.ncol() ; j++) {
      u[i] += XmXinv(i, j) * signs[j] ;
    }
    u[i] *= lambda ;
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
  stepCoef = stepCoef / std::pow(std::max(1, iter + 1 - delay), stepRate) ;
  for(int i = 0 ; i < gradient.length() ; i++) {
    gradient[i] = Xy[i];
    for(int j = 0 ; j < XmX.ncol() ; j++) {
      gradient[i] -= XmX(i, j) * samp[j] ;
    }
    gradient[i] *= stepCoef ;
    if(std::abs(gradient[i]) > gradientBound) {
      gradient[i] = gradientBound * sign_cpp(gradient[i]) ;
    }
  }
}

// [[Rcpp::export]]
double test_computeConditionalMean(NumericVector mu, NumericVector samp,
                                   NumericMatrix XmX, double yvar, int index) {
  return computeConditionalMean_cpp(mu, samp, XmX, yvar, index);
}

// [[Rcpp::export]]
int test_innerZeroCondition(NumericVector samp, NumericVector l,
                            NumericVector u) {
  return innerZeroCondition_cpp(samp, l, u);
}

// [[Rcpp::export]]
int test_innerCheckOneCondition(NumericVector samp, NumericVector thresh) {
  return innerCheckOneCondition_cpp(samp, thresh);
}

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

// [[Rcpp::export]]
NumericVector test_computeOneThreshold(NumericVector signs, double lambda,
                                       NumericMatrix XmXinv) {
  NumericVector u(signs.length());
  computeOneThreshold_cpp(signs, lambda, XmXinv, u);
  return u;
}

// [[Rcpp::export]]
double test_computeDiffThreshold(NumericVector signs, double lambda,
                                 NumericMatrix XmXinv, int coordinate) {
  return computeDiffThreshold_cpp(signs, lambda, XmXinv, coordinate);
}

// [[Rcpp::export]]
double test_sign(double x) {
  return sign_cpp(x);
}

// [[Rcpp::export]]
NumericVector test_boundBeta(NumericVector estimate, NumericVector naive) {
  NumericVector est = clone(estimate);
  boundBeta_cpp(est, naive);
  return est;
}

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

