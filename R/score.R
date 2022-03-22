source('R/collapse.R')

# Computes the continuous ranked probability score (CRPS)
crps <- function(cdf, x, lower = -Inf, upper = Inf) {
  -integrate(function(y, x) (cdf(y) - 1*(y >= x))^2, lower, upper, x = x)$value
}