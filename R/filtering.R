source('R/collapse.R')

# Filtering algorithm:
# Given parameters of the SLDS `params`, and observation sequence `obs`,
# the filtered estimate is computed. The function returns the filtered probability of the
# different states of the hidden discrete variable S, the parameters of the hidden
# continuous variable X and the conditional log-likelihood log(P(y^(t+1) | y^(1:t))) for
# each time step.
filtering <- function(params, obs, method = c('exact', 'gpb1', 'gpb2', 'imm'),
                      display_steps = FALSE) {
  method <- tolower(method[1])
  if (is.null(dim(obs))) obs <- t(obs)  # Turn into row vector
  list2env(params, environment())    # Unpack list: P, A, b, Q, C, d, R, q, nu, Gamma
  theta_H <- mget(c('P', 'A', 'b', 'Q'))  # Model parameters for hidden variables
  theta_O <- mget(c('C', 'd', 'R'))       # Model parameters for observed variable
  M <- nrow(P)             # Number of categories of S
  N <- nrow(A)             # Vector length of X
  t_end <- ncol(obs)
  
  belief_state <- vector(mode = 'list', length = t_end)
  prior_state <- list(
    p_dt     = q,
    W_dt     = matrix(1, M, 1),
    mu_dt    = array(nu, c(N, M, 1)),
    Sigma_dt = array(Gamma, c(N, N, M, 1))
  )
  m0 <- ifelse(method == 'gpb2', 'exact', method)  # At t=0, GPB2 has a different form
  belief_state[[1]] <- cond_step(theta_O, prior_state, obs[, 1], 0, m0)
  if (display_steps) cat("time steps computed: ", 1, '/', t_end, '\n')
  for (t in 1:(t_end - 1)) {
    prior_belief_state <- pred_step(theta_H, belief_state[[t]], t, method)
    belief_state[[t + 1]] <- cond_step(
      theta_O, prior_belief_state, obs[, t + 1], t, method
    )
    if (display_steps) cat("time steps computed: ", t + 1, '/', t_end, '\n')
  }
  belief_state
}

# The forward-propagation step of the filtering algorithm:
# Here the probabilities P(S^(t+1) | y^(1:t)) and P(X^(t+1) |  S^(t + 1), y^(1:t))
# are computed.
pred_step <- function(params, state, t, method) {
  list2env(params, environment())  # Unpack list: P, A, b, Q
  list2env(state, environment())   # Unpack list: p_t, (W_t), mu_t, Sigma_t
  M <- nrow(P)                # Number of categories of S
  N <- nrow(A)                # Vector length of X
  
  if (method == 'exact') {
    J <- M^t
    Jm <- M^(t - 1)
  }
  else if (method %in% c('gpb1', 'imm')) {
    J <- 1
    Jm <- 1
  }
  else if (method == 'gpb2') {
    J <- M
    Jm <- 1
  }
  
  # Initialize new belief state
  W_dt <- matrix(1, M, J)  # Weights where each row sum to 1
  mu_dt <- array(0, c(N, M, J))        # Predictive list of expectation in Gaussian mix
  Sigma_dt <- array(0, c(N, N, M, J))  # Predictive list of variance in Gaussian mix
  p_dt <- c(t(P) %*% p_t)  # Predictive probability of discrete random variable S
  if (method %in% c('exact', 'imm')) {
    W_tilde_dt <- t(P*p_t) / p_dt
  }
  for (k in 1:J) {
    i <- (k - 1) %/% Jm + 1  # 1, 1, ..., 2, 2, ..., 3, 3, ..., M
    j <- (k - 1) %% Jm + 1   # 1, 2, ..., 1, 2, ..., 1, 2, ..., M^(t - 1)
    if (method == 'exact') {
      W_dt[, k] <- W_t[i, j] * W_tilde_dt[, i]  # w_sk, s = 1, ..., M
    }
    for (s in 1:M) {
      # E_k[X^(t+1) | S^(t+1)=s, y^(1:t)] and Var_k[X^(t+1) | S^(t+1)=s, y^(1:t)]
      if (method == 'imm') {
        pp_s <- m_proj(W_tilde_dt[s, ], mu_t[,, j], Sigma_t[,,, j])
        mu_dt[, s, k] <- A[,, s] %*% pp_s$mu + b[, s]
        Sigma_dt[,, s, k] <- Q[,, s] + A[,, s] %*% pp_s$Sigma %*% t(A[,, s])
      }
      else {
        mu_dt[, s, k] <- A[,, s] %*% mu_t[, i, j] + b[, s]
        Sigma_dt[,, s, k] <- Q[,, s] + A[,, s] %*% Sigma_t[,, i, j] %*% t(A[,, s])
      }
    }
  }
  mget(c('p_dt', 'W_dt', 'mu_dt', 'Sigma_dt'))  # Return belief state
}

# The conditioning step of the filtering algorithm:
# Here the probabilities P(S^(t+1) | y^(1:(t+1))) and P(X^(t+1) |  S^(t + 1), y^(1:(t+1))
# are computed.
cond_step <- function(params, state, y, t, method) {
  list2env(params, environment())  # Unpack list: C, d, R
  list2env(state, environment())   # Unpack list: p_dt, (W_dt), mu_dt, Sigma_dt
  M <- nrow(P)             # Number of categories of S
  N <- nrow(A)             # Vector length of X
  if (method == 'exact') {
    J <- M^t
  }
  else if (method %in% c('gpb1', 'imm')) {
    J <- 1
  }
  else if (method == 'gpb2') {
    J <- M
    mu_proj <- array(0, c(N, M, 1))
    Sigma_proj <- array(0, c(N, N, M, 1))
  }
  
  # Initialize new belief state
  W_t <- matrix(0, M, J)    # Weights where each row sum to 1
  p_t <- numeric(M)         # Posterior probability of discrete random variable S
  mu_t <- array(0, c(N, M, J))        # Posterior list of expectation in Gaussian mix
  Sigma_t <- array(0, c(N, N, M, J))  # Posterior list of variance in Gaussian mix
  
  for (s in 1:M) {
    xi_s <- numeric(J)
    for (k in 1:J) {
      K_sk <- t(solve(
        R[,, s] + C[,, s] %*% Sigma_dt[,, s, k] %*% t(matrix(C[,, s], nrow(C))),
        C[,, s] %*% Sigma_dt[,, s, k]
      ))
      mu_t[, s, k] <- mu_dt[, s, k] + K_sk %*% (y - (C[,, s] %*% mu_dt[, s, k] + d[, s]))
      Sigma_t[,, s, k] <- (diag(N) - K_sk %*% C[,, s]) %*% Sigma_dt[,, s, k]
      xi_s[k] <- mvtnorm::dmvnorm(y,
        mean = C[,, s] %*% mu_dt[, s, k] + d[, s],
        sigma = R[,, s] + C[,, s] %*% Sigma_dt[,, s, k] %*% t(matrix(C[,, s], nrow(C))),
        log = FALSE,
        checkSymmetry = TRUE
      )
    }
    p_t[s] <- p_dt[s]*sum(W_dt[s, ] * xi_s)  # Unnormalized
    W_t[s, ] <- (p_dt[s] / p_t[s]) * W_dt[s, ]*xi_s
    if (method == 'gpb2') {
      pp_s <- m_proj(W_t[s, ], mu_t[, s, ], Sigma_t[,, s, ])
      mu_proj[, s, 1] <- pp_s$mu
      Sigma_proj[,, s, 1] <- pp_s$Sigma
    }
  }
  L_t <- sum(p_t)  # Conditional likelihood P(y^(t + 1) | y^(1:t))
  l_t <- log(L_t)  # Log-likelihood
  p_t <- p_t / L_t  # Normalize
  if (method == 'gpb1') {
    # Extract moment-projected µ and Σ (From minimizing KL-divergence)
    pp <- m_proj(p_t, mu_t[,, 1], Sigma_t[,,, 1])  # Projected parameters
    mu_t <- array(pp$mu, c(N, 1, 1))
    Sigma_t <- array(pp$Sigma, c(N, N, 1, 1))
  }
  else if (method == 'gpb2') {
    W_t <- matrix(1, M, 1)
    mu_t <- mu_proj
    Sigma_t <- Sigma_proj
  }
  mget(c('p_t', 'W_t', 'mu_t', 'Sigma_t', 'l_t'))  # Return belief state
}