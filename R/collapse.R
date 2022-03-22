# Moment-projection:
# Here, the first two moments of the Gaussian mixture are matched with a
# non-mixed Gaussian distribution. That is, the Gaussian mixture is Moment-projected
# (M-projected) onto a Gaussian distribution.
m_proj <- function(w, mu_list, Sigma_list) {
  if (is.null(dim(mu_list))) {  # Is univariate
    M <- length(mu_list)  # Number of S categories
    mu_list <- t(mu_list)
    dim(Sigma_list) <- c(1, 1, M)
  }
  N <- nrow(Sigma_list)
  mu <- mu_list %*% w
  Sigma <- matrix(0, N, N)
  for (i in 1:length(w)) {
    Sigma <- Sigma + w[i]*(Sigma_list[,, i] + tcrossprod(mu_list[, i] - mu))
  }
  list(mu = mu, Sigma = Sigma)
}

marginal_x <- function(bs, method) {
  t_end <- length(bs)
  N <- nrow(bs[[1]]$mu_t)
  mu <- matrix(0, N, t_end)
  Sigma <- array(0, c(N, N, t_end))
  if (tolower(method) == 'gpb1') {
    mu[, ] <- vapply(bs, function(x) x$mu_t, numeric(N))
    Sigma[,, ] <- vapply(bs, function(x) x$Sigma_t, numeric(N^2))
  }
  else {
    for (t in 1:t_end) {
      # Unpack list of flatten mix: w, mu_list, Sigma_list
      list2env(flatten_mix(bs[[t]]), environment())
      
      pp <- m_proj(w, mu_list, Sigma_list)
      mu[, t] <- pp$mu
      Sigma[,, t] <- pp$Sigma
    }
  }
  mget(c('mu', 'Sigma'))
}

flatten_mix <- function(bs) {
  w <- c(t(bs$W_t * bs$p_t))
  
  mu_list <- aperm(bs$mu_t, c(1, 3, 2))
  d <- dim(mu_list)
  dim(mu_list) <- c(d[1], d[2]*d[3])
  
  Sigma_list <- aperm(bs$Sigma_t, c(1, 2, 4, 3))
  d <- dim(Sigma_list)
  dim(Sigma_list) <- c(d[1:2], d[3]*d[4])
  
  mget(c('w', 'mu_list', 'Sigma_list'))
}