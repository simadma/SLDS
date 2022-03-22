source('R/utils.R')
library(ggplot2)

M <- 4  # 1=stationary, 2=speeding up, 3=constant velocity, 4=slowing down
N <- 3
L <- 1

q <- c(0.97, 0.01, 0.01, 0.01)

nu_a <- c(0, 1, 0, -1)
sig_a <- c(0, 0.01, 0.01, 0.01)
a_a <- c(0, 0.9, 0.9, 0.9)
b_a <- (1 - a_a)*nu_a
gam_a <- sig_a/sqrt(1 - a_a^2)
dt <- 0.5

nu <- matrix(c(c(0, 0, nu_a[1]),
               c(0, 0, nu_a[2]),
               c(0, 0, nu_a[3]),
               c(0, 0, nu_a[4])),
             N, M)

Gamma <- array(c(c(rep(0, 8), gam_a[1]),
                 c(rep(0, 8), gam_a[2]),
                 c(rep(0, 8), gam_a[3]),
                 c(rep(0, 8), gam_a[4])),
               c(N, N, M))
P <- 0.95*diag(M)
P[cbind(1:M, 1:M %% M + 1)] <- 1-P[[1,1]]
P
A <- array(c(c(1, 0, 0,  0, 0, 0,        0,  0, a_a[1]),
             c(1, 0, 0, dt, 1, 0, 0.5*dt^2, dt, a_a[2]),
             c(1, 0, 0, dt, 1, 0, 0.5*dt^2, dt, a_a[3]),
             c(1, 0, 0, dt, 1, 0, 0.5*dt^2, dt, a_a[4])),
           c(N, N, M))

b <- matrix(c(c(0, 0, b_a[1]),
              c(0, 0, b_a[2]),
              c(0, 0, b_a[3]),
              c(0, 0, b_a[4])),
            N, M)
Q <- array(c(c(rep(0, 8), sig_a[1]),
             c(rep(0, 8), sig_a[2]),
             c(rep(0, 8), sig_a[3]),
             c(rep(0, 8), sig_a[4])),
           c(N, N, M))

C <- array(rep(c(1, 0, 0), M), c(L, N, M))
d <- matrix(0, L, M)
R <- 5^2*array(1, c(L, L, M))
params <- list(P=P, A=A, b=b, Q=Q, C=C, d=d, R=R, q=q, nu=nu, Gamma=Gamma)


## One simulation run
t_end <- 200
time <- dt*(1:t_end)
set.seed(411)
sim <- simulate_slds(params, t_end)
simS_scale <- (sim$S - min(sim$S)) / (max(sim$S) - min(sim$S))*max(sim$X[1,])

## PLOT SIMULATION
par(mfrow = c(1, 1))
s_col <- c(obs=alpha(rgb(0, 0.3, 0.8), 0.5),
           pos=alpha(rgb(1, 0, 0), 0.3),
           vel=alpha(rgb(1, 0, 1), 0.3),
           acc=alpha(rgb(0, 1, 1), 0.3),
           s =alpha(rgb(0, 0, 0), 0.5))
plot(time, sim$Y[1,], type = 'l', col = s_col['obs'],
     ylim = c(0, max(sim$X, sim$Y)),
     main = 'Simulation of SLDS', xlab = 'time', ylab = 'value')
lines(time, sim$X[1,], lty = 1, col = s_col['pos'], lwd = 2)
lines(time, sim$X[2,], lty = 1, col = s_col['vel'], lwd = 2)
lines(time, sim$X[3,], lty = 1, col = s_col['acc'], lwd = 2)
lines(time, simS_scale, type = 'p', col = s_col['s'])
legend('topleft', legend = names(s_col),
       col = s_col, ncol = 3,
       lty = c(1,1,1,1,NA), pch = c(NA,NA,NA,NA,1))


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
  s_bs[[method]] <- vapply(bs[[method]], function(x) x$p_t, numeric(M))
  x_bs[[method]] <- marginal_x(bs[[method]], method)
  l_bs[[method]] <- vapply(bs[[method]], function(x) x$l_t, numeric(1))
}

# PLOT FILTERED ESTIMATES
par(mfrow = c(1, 1))
my_colors <- 1:length(methods)
plot(time, s_bs[[methods[1]]][, ], type = 'l', ylim = c(0, 1), col = my_colors[1],
     main = 'Filtered posterior of S', xlab = 'time', ylab = 'P(S = 2)')
for (i in 2:length(methods)) {
  lines(time, s_bs[[methods[i]]][2, ], type = 'l', col = my_colors[i])
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
