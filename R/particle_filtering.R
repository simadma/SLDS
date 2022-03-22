library(purrr)
library(mvtnorm)

particle_filtering <- function(params, obs, B = 100) {
  if (setequal(names(params), c('P', 'A', 'b', 'Q', 'C', 'd', 'R', 'q', 'nu', 'Gamma'))) {
    list2env(params, environment())    # Unpack list
  } else stop("`params`-list must contain parameters P, A, b, Q, C, d, R, q, nu, Gamma")
  
  M <- nrow(P)             # Number of categories of S
  t_end <- ncol(obs)
  obs_prob <- function(s, x, y) {
    dmvnorm(y, mean = C[,, s] %*% x + d[, s], sigma = as.matrix(R[,, s]))
  }
  sample_s <- function(s_prev) {
    sample.int(M, size = 1, prob = P[s_prev, ])
  }
  sample_x <- function(s, x_prev) {
    c(rmvnorm(1, mean = A[,, s] %*% x_prev + b[, s], sigma = Q[,, s]))
  }
  
  D <- vector('list', length = t_end)  # Samples
  
  # First step
  ss <- sample.int(M, size = B, replace = TRUE, prob = q)
  xs <- map(ss, function(s) c(rmvnorm(1, nu[ , s], Gamma[ , , s])))
  ws <- map2_dbl(ss, xs, obs_prob, y = obs[, 1])
  D[[1]] <- list(ss = ss, xs = xs)
  # The rest
  for (t in seq_len(t_end - 1)) {
    boot <- sample.int(B, replace = TRUE, prob = ws)  # `ws` need not sum to one
    ss <- map_int(D[[t]]$ss[boot], sample_s)
    xs <- map2(ss, D[[t]]$xs[boot], sample_x)
    ws <- map2_dbl(ss, xs, obs_prob, y = obs[, t + 1])
    D[[t + 1]] <- list(ss = ss, xs = xs)
  }
  D
}