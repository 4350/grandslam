library(mvtnorm)

rm(list = ls())
sigma <- matrix(0.50, ncol = 3, nrow = 3)
diag(sigma) <- 1

returns <- rmvnorm(10000, mean = c(0.10, 0.25, 0.00), sigma = sigma)

weights <- c(0.10, 0.30, 0.60)
portfolio_returns <- returns %*% weights

q <- 0.05

# Compute VaR for individual returns
var_limit <- apply(returns, 2, function(r) quantile(r, q))

es_returns <- lapply(seq_along(var_limit),
             function(i) mean(returns[returns[, i] < var_limit[i], i]))

es_returns <- unlist(es_returns)

es_ub <- -sum(weights * es_returns)
es_lb <- -quantile(portfolio_returns, probs = q)

# Portfolio ES
es <- -mean(portfolio_returns[portfolio_returns < -es_lb])

CBD <- (es_ub - es) / (es_ub - es_lb)