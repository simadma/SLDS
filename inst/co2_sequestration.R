source('R/simulate_slds.R')
source('R/filtering.R')
source('R/collapse.R')
source('R/score.R')
library(ggplot2)

M <- 2
N <- 2
L <- 2

q <- c(0.99, 0.01)
nu <- matrix(c(120, 50, 70, 70), N, M)

rho <- 0.8
z <- 20
Gamma <- array(c(c(nu[1, 1]^2,          0,         0,         nu[2, 1]^2),
                 c(nu[1, 2]^2, rep(nu[1, 2]*nu[2, 2]*rho, 2), nu[2, 2]^2))/z^2,
               c(N, N, M))
Gamma
P <- matrix(c(0.99, 0, 0.01, 1), M)
A <- array(rep(c(0.9, 0, 0, 0.9), 2), c(N, N, M))
b <- matrix(0, N, M)
Q <- array(0, c(N, N, M))
for (s in 1:M) {
  b[, s] <- (diag(N) - A[,, s]) %*% nu[, s]
  Q[,, s] <- Gamma[,, s] - A[,, s] %*% Gamma[,, s] %*% t(A[,, s])
}

C <- array(rep(c(1, 0, 0, 1), 2), c(L, N, M))
d <- matrix(0, L, M)
R <- 15^2*array(rep(c(1, 0, 0,  1), 2), c(L, L, M))
params <- list(P=P, A=A, b=b, Q=Q, C=C, d=d, R=R, q=q, nu=nu, Gamma=Gamma)


## One simulation run
t_end <- 100
set.seed(42)
sim <- simulate_slds(params, t_end)

## PLOT SIMULATION
par(mfrow = c(1, 1))
s_col <- c(y1=alpha(rgb(0, 0.3, 0.8), 0.5),
           y2=alpha(rgb(0, 0.8, 0.3), 0.5),
           x1=alpha(rgb(1, 0, 0), 0.3),
           x2=alpha(rgb(1, 0, 1), 0.3),
           s =alpha(rgb(0, 0, 0), 0.5))
plot(sim$Y[1,], type = 'l', col = s_col['y1'],
     ylim = c(0, max(sim$X, sim$Y)),
     main = 'Simulation of SLDS', xlab = 'time', ylab = 'value')
lines(sim$X[1,], lty = 1, col = s_col['x1'], lwd = 2)
lines(sim$Y[2,], lty = 1, col = s_col['y2'])
lines(sim$X[2,], lty = 1, col = s_col['x2'], lwd = 2)
lines(ifelse(sim$S == 1, 0, max(sim$Y, sim$X)), type = 'p', col = s_col['s'])
legend('bottomright', legend = names(s_col),
       col = s_col, ncol = 3,
       lty = c(1,1,1,1,NA), pch = c(NA,NA,NA,NA,1))

# saveRDS(params, file = "./data/co2_params.rds")
# saveRDS(sim, file = "./data/co2_sim.rds")
## COMPUTE FILTERED ESTIMATES
y <- sim$Y  # Evidence

bs <- list()
s_bs <- list()
x_bs <- list()
l_bs <- list()
if (t_end <= 13) {
  methods <- c('GPB1', 'IMM', 'GPB2', 'Exact')
} else {
  methods <- c('GPB1', 'IMM', 'GPB2')
}

for (method in methods) {
  bs[[method]] <- filtering(params, y, method)
  s_bs[[method]] <- vapply(bs[[method]], function(x) x$p_t, numeric(2))
  x_bs[[method]] <- marginal_x(bs[[method]], method)
  l_bs[[method]] <- vapply(bs[[method]], function(x) x$l_t, numeric(1))
}

# PLOT FILTERED ESTIMATES
par(mfrow = c(1, 1))
my_colors <- 1:length(methods)
plot(s_bs[[methods[1]]][2, ], type = 'l', ylim = c(0, 1), col = my_colors[1],
     main = 'Filtered posterior of S', xlab = 'time', ylab = 'P(S = 2)')
for (i in 2:length(methods)) {
  lines(s_bs[[methods[i]]][2, ], type = 'l', col = my_colors[i])
}
legend('topleft', legend = methods, col = my_colors, lty = 1,
       y.intersp = 0.3)
abline(v = match(2, sim$S), col = 'red', lty = 2)

plot(x_bs[[methods[1]]]$mu[1, ], type = 'l', col = my_colors[1],
     main = 'Filtered posterior of X', xlab = 'time', ylab = 'E[X^t | y^(1:t)]')
for (i in 2:length(methods)) {
  lines(x_bs[[methods[i]]]$mu[1, ], type = 'l', col = my_colors[i])
}
abline(v = match(2, sim$S), col = 'red', lty = 2)
legend('topleft', legend = methods, col = my_colors, lty = 1,
       y.intersp = 0.3)


##########################################################################################
## PLOT CONDITIONAL LOG-LIKELIHOOS
par(mfrow = c(1, 1))
plot(l_bs[[methods[1]]], type = 'l', col = my_colors[1],
     main = expression("Conditional log-likelihood:" ~log~p(y^{(t+1)} ~ '|' ~ y^{(1:t)})),
     ylim = c(min(sapply(l_bs, min)), max(sapply(l_bs, max))),
     xlab = 'time', ylab = 'log-likelihood')
for (i in 2:length(methods)) {
  lines(l_bs[[methods[[i]]]], col = my_colors[i])
}
legend('bottomleft', legend = methods, col = my_colors, lty = 1)
##########################################################################################
