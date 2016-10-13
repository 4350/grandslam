library(devtools)
library(mvtnorm)
load_all('wimbledon')

sigma <- matrix(0.50, 6, 6)
diag(sigma) <- 1

set.seed(403)
x <- rmvnorm(2500, mean = c(0, 0, 0, 0, 0, 0),
             sigma = sigma)

u <- pnorm(x)

dist.params <- list(df = 8, skew = rep(0.10, 6))
tic('Run Model')
results <- wimbledon::dc.run.model(u, dist.params, alpha = 0.06, beta = 0.91)
toc()

results$Omega
head(results$shocks, 2)
wimbledon::dc.ll(dist.params, results)