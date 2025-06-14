% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lasso.R
\name{lassoMLE}
\alias{lassoMLE}
\title{Maximum Likelihood Inference for Models Selected by the LASSO}
\usage{
lassoMLE(y, X, lambda = "lambda.min", ysig = NULL, lassoFit = NULL,
  delay = 50, optimSteps = 1000, sampSteps = 2000, stepCoef = 0.001,
  stepRate = 0.85, method = c("exact", "selected"))
}
\arguments{
\item{y}{response variable.}

\item{X}{input matrix of dimensions nobs x nvars.}

\item{lambda}{regularization parameter value to be used with
\code{\link[glmnet]{glmnet}}.}

\item{ysig}{the residual standard error. If \code{NULL} then the
cross-validated lasso estimate will be used.}

\item{lassoFit}{a lasso fit to the data, of the sort returned by
\code{\link[glmnet]{glmnet}} or \code{\link[glmnet]{cv.glmnet}}. If
NULL, the lasso fit will be computed within the routine.}

\item{delay}{number of optimization steps to take before decreasing
the gradient step size.}

\item{optimSteps}{number of stochastic gradient steps to take.}

\item{sampSteps}{number of samples to take in order to compute the confidence
intervals.}

\item{stepCoef}{fixed step size for stochastic gradient.}

\item{stepRate}{the rate at which to decrease the step size of the
stochastic gradient.}

\item{method}{optimization method, "selected" is faster than "exact"
but may not converge to the exact MLE. See description for details.}

\item{verbose}{logical; if \code{TRUE} a progress bar is displayed during
 optimisation and sampling.}
}
\value{
\code{lassoMLE} returns an object of class \code{lassoMLE} which
  contains the following variables:

  * \code{lassoFit} the \code{\link[glmnet]{glmnet}} or
  \code{\link[glmnet]{cv.glmnet}} objects describing the solution path
  of the lasso.

  * \code{activeSet} the models selected by the lasso.

  * \code{conditionalBeta} the conditional MLE for the model selected
  by the lasso.

  * \code{lassoBeta} the lasso coefficients estimates.

  * \code{lassoysig} the cross-validated lasso residual standard
  error estimate.

  * \code{wald_CI} the conditional post-selection confidence intervals.

  * \code{solution_path} the solution path of the stochastic gradient method.
  This is not the lasso solution path!

  * \code{coef_sample} a matrix of samples from the post-selection distribution
  of the refitted regression coefficients. These were used to compute the conditional
  confidence intervals.
}
\description{
\code{lassoMLE} computes the conditional MLE for models
selected by the lasso using a stochastic gradient/contrastive divergence
approach. Post-selection confidence intervals are computed based on the
theoretical approximate distribution of the conditional MLE.
}
\details{
The routine computes the conditional MLE for models selected
by the lasso as well as post-selection confidence intervals which are
based on the approximate distribution of the conditional MLE. The routine
uses the \code{\link[glmnet]{cv.glmnet}} to perform model selection and then
estimates the regression coefficietns using a stochastic gradient approach.

The distinction between the `exact' and `selected' estimation methods is that
for exact all the variables relevant to the selection of the lasso model are sampled,
including variables related to the set of variables which were not included in the model.
Because the number of unselected variables tends to be far larger than the number
of variables included in the model, the `selected' method which only samples the refitted
regression coefficients tends to be much faster. The two methods tend to produce very similar
results.
}
