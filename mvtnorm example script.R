## Example demonstrating truncNormMLE in a terminal friendly way
suppressPackageStartupMessages({
  library(selectiveMLE)
  library(knitr)
})

ascii_hist <- function(x, bins = 10, width = 50) {
  h <- hist(x, breaks = bins, plot = FALSE)
  max_count <- max(h$counts)
  for (i in seq_along(h$counts)) {
    bar <- if (max_count > 0) {
      paste(rep("#", round(width * h$counts[i] / max_count)), collapse = "")
    } else ""
    cat(sprintf("%8.2f | %s\n", h$mids[i], bar))
  }
}
# Simulation Parameters --------------
set.seed(123)
p <- 100
rho <- 0.3
sigma <- matrix(rho, nrow = p, ncol = p)
diag(sigma) <- 1
threshold <- 1.65
plarge <- 0.2
signalsd <- 2
nullsd <- 0.0000005
cialpha <- 0.05

message("Generating data ...")
ind <- rbinom(p, 1, plarge)
mu <- (1 - ind) * rnorm(p, mean = 0, sd = nullsd) + ind * rnorm(p, mean = 0, sd = signalsd)
mu <- rep(0, p)
mu[1:20] <- rnorm(20, 0, 2)
y <- rep(0, p)
while (all(y < abs(threshold))) {
  y <- as.vector(mvtnorm::rmvnorm(1, mu, sigma))
}
selected <- abs(y) > threshold

message("Computing conditional MLE ...")
fit <- truncNormMLE(y, sigma, threshold, cialpha = cialpha,
                    maxiter = 500)

# Solution path (first selected parameter)
path <- fit$solutionPath
path <- path[, selected, drop = FALSE]
message("Solution path for parameter 1:")
ascii_hist(path[, 1])

# Summarise estimates and confidence intervals
true <- mu[selected]
naive <- y[selected]
lci <- naive + qnorm(cialpha / 2) * sqrt(diag(sigma)[selected])
uci <- naive + qnorm(1 - cialpha / 2) * sqrt(diag(sigma)[selected])
naivecover <- mean(lci < true & uci > true)
conditional <- fit$mle
condcover <- mean(fit$CI[, 1] < true & fit$CI[, 2] > true)

tab <- data.frame(
  estimate = conditional,
  CI_lower = fit$CI[, 1],
  CI_upper = fit$CI[, 2],
  naive_est = naive
)
message("Coverage (naive, conditional):")
print(c(naivecover, condcover))
message("\nEstimates and CIs:")
knitr::kable(head(tab))


