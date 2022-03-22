# Simulation of an SLDS:
# The parameters are
# P, a transition probability matrix for the hidden S variable,

# A, a matrix A_s for each s for the CLG N(A_s x + b_s, Q_s)
# b, a vector b_s for each s for the CLG N(A_s x + b_s, Q_s)
# Q, a matrix Q_s for each s for the CLG N(A_s x + b_s, Q_s)

# C, a matrix C_s for each s for the CLG N(C_s x + d_s, R_s)
# d, a vector d_s for each s for the CLG N(C_s x + d_s, R_s)
# R, a matrix R_s for each s for the CLG N(C_s x + d_s, R_s)

# q, initial probabilities of S
# nu, a vector nu_s for each s for the CLG N(nu_s, Q_s)
simulate_slds <- function(params, n) {
  list2env(params, environment())  # Unpack list: P, A, b, Q, C, d, R, q, nu, Gamma
  M <- nrow(P)
  N <- nrow(A)
  L <- nrow(C)
  
  S <- numeric(n)
  X <- matrix(0, N, n)
  Y <- matrix(0, L, n)
  
  # First step
  S[1] <- sample(M, 1, prob = q)
  X[, 1] <- mvtnorm::rmvnorm(1, nu[, S[1]], as.matrix(Gamma[,, S[1]]))
  Y[, 1] <- mvtnorm::rmvnorm(1, C[,, S[1]] %*% X[, 1] + d[, S[1]], as.matrix(R[,, S[1]]))
  # Iterate
  for (t in 1:(n - 1)) {
    S[t + 1] <- sample(M, 1, prob = P[S[t], ])
    X[, t + 1] <- mvtnorm::rmvnorm(1,
      A[,, S[t + 1]] %*% X[, t] + b[, S[t + 1]],
      as.matrix(Q[,, S[t + 1]])
    )
    Y[, t + 1] <- mvtnorm::rmvnorm(1,
      C[,, S[t + 1]] %*% X[, t + 1] + d[, S[t + 1]],
      as.matrix(R[,, S[t + 1]])
    )
  }
  list(S = S, X = X, Y = Y)
}