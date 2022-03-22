source('R/simulate_slds.R')
source('R/filtering.R')
source('R/collapse.R')
source('R/score.R')
library(ggplot2)

M <- 2
N <- 2
L <- 2

P <- matrix(c(0.99, 0, 0.01, 1), M)
A <- array(rep(c(0.9, 0, 0, 0.9), 2), c(N, N, M))
b <- matrix(c(10,  2, 8, 8), N, M)
rho <- 0.4
Q <- array(c(c(5^2, 0, 0, 1), 4^2*c(1, rho, rho, 1)), c(N, N, M))

C <- array(rep(c(1, 0, 0, 1), 2), c(L, N, M))
d <- matrix(0, L, M)
R <- array(rep(c(20^2, 0, 0,  20^2), 2), c(L, L, M))

q <- c(0.99, 0.01)

nu <- matrix(c(100, 20, 80, 80), N, M)

t_end <- 100
set.seed(42)
sim <- simulate_slds(mget(c('P', 'A', 'b', 'Q', 'C', 'd', 'R', 'q', 'nu')), t_end)
y <- sim$Y

## PLOT SIMULATION
par(mfrow = c(1, 1))
plot(sim$Y[1,], type = 'l', col = alpha(rgb(0, 0.3, 0.8), 0.5),
     ylim = c(min(sim$X, sim$Y), max(sim$X, sim$Y)),
     main = 'Simulation of SLDS', xlab = 'time', ylab = 'value')
lines(sim$X[1,], lty = 1, col = alpha(rgb(1, 0, 0), 0.3), lwd = 2)
lines(sim$Y[2,], lty = 1, col = alpha(rgb(0, 0.8, 0.3), 0.5))
lines(sim$X[2,], lty = 1, col = alpha(rgb(1, 0, 1), 0.3), lwd = 2)
lines(ifelse(sim$S == 1, min(sim$X), max(sim$X)), type = 'p', col = 'gray35')
legend('bottomright', legend = c('Y', 'X', 'S'), col = c('black', 'gray35', 'gray35'),
       lty = c(1,2,NA), pch = c(NA,NA,1))

## COMPUTE FILTERED ESTIMATES
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
  bs[[method]] <- filtering(mget(c('P', 'A', 'b', 'Q', 'C', 'd', 'R', 'q', 'nu')), y, method)
  s_bs[[method]] <- vapply(bs[[method]], function(x) x$p_t, numeric(2))
  x_bs[[method]] <- marginal_x(bs[[method]], method)
  l_bs[[method]] <- vapply(bs[[method]], function(x) x$l_t, numeric(1))
}

# PLOT FILTERED ESTIMATES
par(mfrow = c(2, 1))
my_colors <- 1:length(methods)
plot(s_bs[[methods[1]]][2, ], type = 'l', ylim = c(0, 1), col = my_colors[1],
     main = 'Filtered posterior of S', xlab = 'time', ylab = 'P(S = 2)')
for (i in 2:length(methods)) {
  lines(s_bs[[methods[i]]][2, ], type = 'l', col = my_colors[i])
}
legend('topleft', legend = methods, col = my_colors, lty = 1,
       y.intersp = 0.3)


plot(x_bs[[methods[1]]]$mu[1, ], type = 'l', col = my_colors[1],
     main = 'Filtered posterior of X', xlab = 'time', ylab = 'E[X^t | y^(1:t)]')
for (i in 2:length(methods)) {
  lines(x_bs[[methods[i]]]$mu[1, ], type = 'l', col = my_colors[i])
}
legend('topleft', legend = methods, col = my_colors, lty = 1,
       y.intersp = 0.3)

##########################################################################################
# SCORES

# crps_Gauss <- function(x, mu = 0, sigma = 1) {
#   sigma*(1/sqrt(pi) - 2*dnorm((x - mu)/sigma) - (x - mu)/sigma*(2*pnorm((x - mu)/sigma) - 1))
# }
# 
# source('R/score.R')
# cdf <- cdf_mix(bs[[methods[1]]][[1]], methods[1])()
# cdf(1:3)
# 
# crps <- function(cdf, x, lower = -Inf, upper = Inf) {
#   -integrate(function(y, x) (cdf(y) - 1*(y >= x))^2, lower, upper, x = x)$value
# }
# 
# crps(cdf_mix(bs[[methods[1]]][[1]], methods[1])(), y[1, t])
# cdf_mix(bs[[methods[1]]][[1]], methods[1])()
# 
# 
# vapply(1, function(x) mvtnorm::pmvnorm(-Inf, x, mean = 0, sigma = as.matrix(1)), numeric(1))
# cdf <- cdf_mix(bs[[methods[1]]][[5]], methods[1])()


# crps(cdf_mix(bs[[methods[1]]][[5]], methods[1])(), y[1, t])
# score <- vapply(1:length(bs$GPB1),
#   function(t) crps(
#     cdf_mix(bs[[methods[1]]][[t]], methods[1])(),
#     y[1, t]
#   ),
#   numeric(1)
# )

# bs$GPB1[[1]]
# score <- vapply(x, function(x) crps(function(x) pnorm(x, 2, 2), x), numeric(1))
# score2 <- crps_Gauss(x,2,2)
# 
# plot(x, score, type = 'l')
# lines(x, score2, lty = 2, col='red')
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

## COMPUTE RMSE
set.seed(131)
nsim <- 100
sim_bs <- list()
s_rmse <- list()
x_rmse <- list()
for (method in methods) {
  sim_bs[[method]] <- simulate_bs(bs[[method]], method, nsim)
  
  s_se <- (t(sim_bs[[method]]$s) - sim$S)^2
  s_rmse[[method]] <- sqrt(rowMeans(s_se))
  
  x_se <- aperm((aperm(sim_bs[[method]]$x, c(1, 3, 2)) - rep(sim$X, nsim))^2, c(1, 3, 2))
  x_rmse[[method]] <- sqrt(apply(x_se, c(1, 3), mean))
}

## PLOT RMSE
par(mfrow=c(2,1))
plot(s_rmse[[methods[1]]], type = 'l', col = my_colors[1],
     ylim = c(0, max(sapply(s_rmse, max))),
     main = paste0("RMSE of S from n=", nsim, " simulations"), xlab = 'time', ylab = 'RMSE')
for (i in 1:length(methods)) {
  lines(s_rmse[[methods[i]]], col = my_colors[i])
}
legend('topleft', legend = methods, col = my_colors, lty = 1,
       y.intersp = 0.3)

plot(x_rmse[[methods[1]]][1, ], type = 'l', col = my_colors[1],
     ylim = c(0, max(sapply(x_rmse, max))),
     main = paste0("RMSE of X from n=", nsim, " simulations"), xlab = 'time', ylab = 'RMSE')
for (i in 1:length(methods)) {
  lines(x_rmse[[methods[i]]][1, ], col = my_colors[i])
}
legend('topleft', legend = methods, col = my_colors, lty = 1,
       y.intersp = 0.3)


##########################################################################################

s_means <- colMeans(sim_bs$GPB1$s)
x_means <- matrix(vapply(1:t_end, function(t) rowMeans(matrix(sim_bs$GPB1$x[,, t], N)), numeric(N)), N)
x_sds <- matrix(vapply(1:t_end, function(t) sqrt(rowMeans((matrix(sim_bs$GPB1$x[,, t], N) - x_means[, t])^2)), numeric(N)), N)


par(mfrow = c(2, 1))
plot(s_means - 1, type = 'l')
lines(s_bs$GPB1[2, ], col = 'red')


plot(x_means[1, ], type='l')
lines(x_means[1, ] - 1.96*x_sds[1, ], lty = 2)
lines(x_means[1, ] + 1.96*x_sds[1, ], lty = 2)

lines(x_bs$GPB1$mu[1, ], col = 'red')
lines(x_bs$GPB1$mu[1, ] - 1.96*sqrt(x_bs$GPB1$Sigma[1, 1, ]), col = 'red', lty = 2)
lines(x_bs$GPB1$mu[1, ] + 1.96*sqrt(x_bs$GPB1$Sigma[1, 1, ]), col = 'red', lty = 2)
lines(sim$X[1, ], col = 'blue')
lines(sim$Y[1, ], col = 'cyan', lty = 2)



apply(s_bs$exact, 2, function(x) which.max(x))
apply(s_bs$imm, 2, function(x) which.max(x))

