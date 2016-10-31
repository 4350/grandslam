library(foreach)
library(tictoc)

rm(list = ls())
NAME <- 'full_indep'
gamma = 5

load(sprintf('data/derived/distributions/%s_10000.RData', NAME))
distribution <- exp(distribution) - 1

load('data/derived/weekly-full.RData')

e_crra_utility <- function(w, gamma, distribution) {
  portfolio <- distribution %*% w
  mean(((1 + portfolio) ^ (1 - gamma)) / (1 - gamma))
}

optimize_utility_t <- function(gamma, distribution_t) {
  N <- ncol(distribution_t)
  
  # Need to start in interior
  theta <- rep(1 / N - 0.01, N)
  
  ui <- rbind(
    diag(1, N),        # lower bound
    rep(-1, N)         # sum to one
  )
  ci <- c(rep(0, N), -1)
  
  constrOptim(
    theta,
    e_crra_utility,
    grad = NULL,
    ui = ui,
    ci = ci,
    
    control = list(
      fnscale = -1
    ),
    
    gamma = gamma,
    distribution = distribution_t
  )
}

times <- 1:dim(distribution)[3]

# Pick the out of sample returns, and then the subset we're testing here
realized <- tail(df[, -1], nrow(distribution))
realized <- head(realized, length(times))

results <- foreach (t = times) %do% {
  tic(sprintf('CRRA (gamma = %d; t = %d', gamma, t))
  o <- optimize_utility_t(gamma, distribution[,, t])
  toc()
  
  list(
    weights = o$par,
    crra = o$value,
    realized = sum(realized[t, ] * o$par)
  )
}

save(results, file = sprintf('data/derived/crra/simple_%d_%s.RData', gamma, NAME))