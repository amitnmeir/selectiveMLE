# selectiveMLE

`selectiveMLE` is an R package that implements conditional maximum likelihood
estimation for models that have been fit after a selection step.  It
currently provides tools for

* `lassoMLE()` – inference for linear models selected by the LASSO; and
* `truncNormMLE()` – inference for selected coordinates of a multivariate
  normal vector.

The package relies on `glmnet` and a set of C++ routines (via `Rcpp`) to
perform stochastic gradient based optimisation.  It is intended for research
purposes and does not currently appear on CRAN.

## Installation

Clone this repository and build the package with `R CMD INSTALL` or
from within R using `devtools`:

```r
# install.packages("devtools")
library(devtools)
install_local("path/to/selectiveMLE")
```

Alternatively, use the provided `Makefile`:

```bash
make install
```

This calls `R CMD INSTALL` after running `make clean`, which removes any
compiled artifacts from previous builds. You can invoke `make clean` on its
own to clear these files without installing the package.

## Quick example

```r
library(selectiveMLE)
set.seed(123)
X <- matrix(rnorm(100 * 10), 100, 10)
beta <- c(2, -1.5, rep(0, 8))
y <- as.numeric(X %*% beta + rnorm(100))

# compute conditional MLE for the model chosen by the lasso
fit <- lassoMLE(y, X, lambda = "lambda.min", verbose = FALSE)
fit$conditionalBeta
fit$wald_CI
```

For truncated normal means:

```r
sigma <- diag(5)
y <- rnorm(5, mean = 0, sd = 1.5)
res <- truncNormMLE(y, sigma, threshold = 1, maxiter = 100, verbose = FALSE)
res$mle
res$CI
```

More elaborate demonstrations can be found in the scripts
`lassoMLE example script.R` and `mvtnorm example script.R` at the top of
this repository. The latter now logs progress and prints results with
`knitr::kable` and ASCII histograms produced by `txtplot` so that it can be
run easily in a terminal session.

## Running tests

The package uses the `testthat` framework.  From the package directory run:

```r
library(testthat)
test_dir("tests/testthat")
```

This will execute the basic unit tests for both exported functions.

