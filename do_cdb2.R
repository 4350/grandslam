# Setup ------------------------------------------------------------------

rm(list = ls())

library(tictoc)
library(devtools)
load_all('australian')
load('data/derived/weekly-estim.RData')

MODEL_NAME <- 'full_dynamic_std_10000'

optim_cdb_constrOptim <- function(weights, fn) {
  N <- length(weights)
  
  # Don't optimize on the first weight, as it is restricted by the
  # sum(weights) == 1 condition.
  weights1 <- weights[-1]
  
  # Call constrOptim with different methods
  optimize <- function(weights1, method = 'Nelder-Mead') {
    constrOptim(
      weights1,
      function(w) -fn(c(1 - sum(w), w)),
      grad = NULL,
      
      ui = rbind(
        diag(1, N - 1),   # All greater than zero
        rep(-1, N - 1)    # Sum less than one
      ),
      
      ci = c(rep(0, N - 1), -1),
      
      control = list(
        maxit = 50000
      ),
      
      method = method
    )
  }
  
  # First, try Nelder-Mead optimization
  op <- optimize(weights1)
  
  # If that fails, try SANN
  if (op$convergence > 0) {
    op <- optimize(op$par, method = 'SANN')
  }
  
  list(
    weights = c(1 - sum(op$par), op$par),
    cdb = -op$value
  )
}

do_best_cdb <- function(model_name, strategy, selectors, q = 0.05) {
  load(sprintf('data/derived/distributions/%s.RData', model_name))
  
  # Restrict to the given subset
  colnames(distribution) <- c('Mkt.RF', 'HML', 'SMB', 'Mom', 'RMW', 'CMA')
  distribution <- distribution[, selectors, ]
  
  times <- 1:dim(distribution)[3]
  N <- ncol(distribution)
  
  # Output: Optimal Weights and CDB at optimal weights
  weights <- matrix(NA, ncol = N, nrow = length(times))
  colnames(weights) <- selectors
  cdb <- rep(NA, length(times))
  
  for (t in times) {
    fn <- cdb_fn(q, distribution[,, t])
    
    tic(sprintf("Optimal CDB for t = %d", t))
    op <- optim_cdb_constrOptim(fn = fn, weights = rep(1 / N, N))
    toc()
    
    weights[t, ] <- op$weights
    cdb[t] <- op$cdb
  }
  
  # Get dates, realized returns, cdb
  metrics <- portfolio_metrics(weights, distribution, q, selectors)
  
  # Make results list
  
  results <- c(
    metrics,
    list(
      weights = weights,
      cdb_check = cdb
    )
  )
  
  file <- sprintf('data/derived/cdb/constrOptim_q%.0f_%s_%s.RData', 100 * q, model_name, strategy)
  save(results, file = file)
}

do_ew_cdb <- function(model_name, strategy, selectors, q = 0.05, df.realized = df.estim) {
  load(sprintf('data/derived/distributions/%s.RData', model_name))
  
  # Restrict to universe
  colnames(distribution) <- c('Mkt.RF', 'HML', 'SMB', 'Mom', 'RMW', 'CMA')
  distribution <- distribution[, selectors, ]
  
  times <- 1:dim(distribution)[3]
  N <- ncol(distribution)
  
  # Output: Optimal Weights and CDB at EW
  weights <- matrix(NA, ncol = N, nrow = length(times))
  colnames(weights) <- selectors
  cdb <- rep(NA, length(times))
  
  for (t in times) {
    fn <- cdb_fn(q, distribution[,, t])
    
    weights[t, ] <- rep(1 / N, N)
    cdb[t] <- fn(weights[t, ])
  }
  
  # Get dates, realized returns, cdb
  metrics <- portfolio_metrics(weights, distribution, q, selectors)
  
  # Make results list
  
  results <- c(
    metrics,
    list(
      weights = weights,
      cdb_check = cdb
    )
  )
  
  file <- sprintf('data/derived/cdb/ew_q%.0f_%s_%s.RData', 100 * q, model_name, strategy)
  save(results, file = file)
}

# CDB Optimization -------------------------------------------------------

do_best_cdb(MODEL_NAME, '5F',          c('Mkt.RF', 'HML', 'SMB', 'RMW', 'CMA'))
do_best_cdb(MODEL_NAME, '5F_EXCL_HML', c('Mkt.RF',        'SMB', 'RMW', 'CMA'))
do_best_cdb(MODEL_NAME, '5F_EXCL_CMA', c('Mkt.RF', 'HML', 'SMB', 'RMW'       ))
do_best_cdb(MODEL_NAME, '5F_EXCL_RMW', c('Mkt.RF', 'HML', 'SMB',        'CMA'))

do_best_cdb(MODEL_NAME, '6F',          c('Mkt.RF', 'HML', 'SMB', 'Mom', 'RMW', 'CMA'))
do_best_cdb(MODEL_NAME, '6F_EXCL_HML', c('Mkt.RF',        'SMB', 'Mom', 'RMW', 'CMA'))
do_best_cdb(MODEL_NAME, '6F_EXCL_CMA', c('Mkt.RF', 'HML', 'SMB', 'Mom', 'RMW'       ))
do_best_cdb(MODEL_NAME, '6F_EXCL_RMW', c('Mkt.RF', 'HML', 'SMB', 'Mom',        'CMA'))

# EW CDB -----------------------------------------------------------------

do_ew_cdb(MODEL_NAME, '5F',          c('Mkt.RF', 'HML', 'SMB', 'RMW', 'CMA'))
do_ew_cdb(MODEL_NAME, '5F_EXCL_HML', c('Mkt.RF',        'SMB', 'RMW', 'CMA'))
do_ew_cdb(MODEL_NAME, '5F_EXCL_CMA', c('Mkt.RF', 'HML', 'SMB', 'RMW'       ))
do_ew_cdb(MODEL_NAME, '5F_EXCL_RMW', c('Mkt.RF', 'HML', 'SMB',        'CMA'))

do_ew_cdb(MODEL_NAME, '6F',          c('Mkt.RF', 'HML', 'SMB', 'Mom', 'RMW', 'CMA'))
do_ew_cdb(MODEL_NAME, '6F_EXCL_HML', c('Mkt.RF',        'SMB', 'Mom', 'RMW', 'CMA'))
do_ew_cdb(MODEL_NAME, '6F_EXCL_CMA', c('Mkt.RF', 'HML', 'SMB', 'Mom', 'RMW'       ))
do_ew_cdb(MODEL_NAME, '6F_EXCL_RMW', c('Mkt.RF', 'HML', 'SMB', 'Mom',        'CMA'))