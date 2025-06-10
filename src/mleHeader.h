#include <Rcpp.h>
using namespace Rcpp;

/*
 * Simple header exposing helper functions used across the C++ source
 * files.  These functions implement truncated normal sampling and
 * conditional means required for the selective MLE algorithms.
 */
/// Sample from a univariate normal distribution truncated below at `threshold`.
double sampleUnivTruncNorm(double mu, double sd, double threshold);

/// Compute the conditional mean of coordinate `index` given the others.
double computeConditionalMean(NumericVector mu,
                              NumericVector samp,
                              const NumericMatrix XmX,
                              double yvar,
                              int index);
