# Setup ------------------------------------------------------------------

rm(list = ls())

library(tictoc)
library(alabama)
library(devtools)
load_all('australian')

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

# optim_cdb_constrOptim.nl <- function(weights, fn) {
#   op <- constrOptim.nl(
#     # Start with EW portfolio; don't optimize on the "first" parameter
#     rep(1 / N, N - 1),
#     function(w) fn(c(1 - sum(w), w)),
#     
#     # Greater than zero; sum le 1
#     hin = function(w) c(w, 1 - sum(w)),
#     
#     # Don't talk so much!!!
#     control.outer = list(
#       trace = FALSE
#     )
#   )
# }

do_best_cdb <- function(model_name, strategy, selectors, q = 0.05) {
  load(sprintf('data/derived/distributions/%s.RData', model_name))
  distribution_simple <- distribution
  rm(distribution)
  
  # Restrict to the given subset
  colnames(distribution_simple) <- c('Mkt.RF', 'HML', 'SMB', 'Mom', 'RMW', 'CMA')
  distribution_simple <- distribution_simple[, selectors, ]
  
  times <- 1:dim(distribution_simple)[3]
  N <- ncol(distribution_simple)
  
  # Output: Optimal Weights and CDB at optimal weights
  weights <- matrix(NA, ncol = N, nrow = length(times))
  colnames(weights) <- selectors
  cdb <- rep(NA, length(times))
  
  for (t in times) {
    fn <- cdb_fn(q, distribution_simple[,, t])
    
    tic(sprintf("Optimal CDB for t = %d", t))
    op <- optim_cdb_constrOptim(fn = fn, weights = rep(1 / N, N))
    toc()
    
    weights[t, ] <- op$weights
    cdb[t] <- op$cdb
  }
  
  cdb_results <- list(
    cdb = cdb,
    weights = weights
  )
  
  file <- sprintf('data/derived/cdb/constrOptim_q%.0f_%s_%s.RData', 100 * q, model_name, strategy)
  save(cdb_results, file = file)
}

do_ew_cdb <- function(model_name, strategy, selectors, q = 0.05) {
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
  
  cdb_results <- list(
    cdb = cdb,
    weights = weights
  )
  
  file <- sprintf('data/derived/cdb/ew_q%.0f_%s_%s.RData', 100 * q, model_name, strategy)
  save(cdb_results, file = file)
}

# CDB Optimization -------------------------------------------------------

MODEL_NAME <- 'full_dynamic_std_10000'

do_best_cdb(MODEL_NAME, '5F',          c('Mkt.RF', 'HML', 'SMB', 'RMW', 'CMA'))
do_best_cdb(MODEL_NAME, '5F_EXCL_HML', c('Mkt.RF',        'SMB', 'RMW', 'CMA'))
do_best_cdb(MODEL_NAME, '5F_EXCL_CMA', c('Mkt.RF', 'HML', 'SMB', 'RMW'       ))
do_best_cdb(MODEL_NAME, '5F_EXCL_RMW', c('Mkt.RF', 'HML', 'SMB',        'CMA'))

do_best_cdb(MODEL_NAME, '6F',          c('Mkt.RF', 'HML', 'SMB', 'Mom', 'RMW', 'CMA'))
do_best_cdb(MODEL_NAME, '6F_EXCL_HML', c('Mkt.RF',        'SMB', 'Mom', 'RMW', 'CMA'))
do_best_cdb(MODEL_NAME, '6F_EXCL_CMA', c('Mkt.RF', 'HML', 'SMB', 'Mom', 'RMW'       ))
do_best_cdb(MODEL_NAME, '6F_EXCL_RMW', c('Mkt.RF', 'HML', 'SMB', 'Mom',        'CMA'))

do_best_cdb(MODEL_NAME, '3F_HML',     c('Mkt.RF', 'HML', 'SMB'))
do_best_cdb(MODEL_NAME, '3F_HML_RMA', c('Mkt.RF', 'HML', 'SMB', 'RMW'))
do_best_cdb(MODEL_NAME, '3F_CMA',     c('Mkt.RF', 'CMA', 'SMB'))
do_best_cdb(MODEL_NAME, '3F_CMA_RMW', c('Mkt.RF', 'CMA', 'SMB', 'RMW'))

do_best_cdb(MODEL_NAME, '4F_HML',     c('Mkt.RF', 'HML', 'SMB', 'Mom'))
do_best_cdb(MODEL_NAME, '4F_HML_RMW', c('Mkt.RF', 'HML', 'SMB', 'Mom', 'RMW'))
do_best_cdb(MODEL_NAME, '4F_CMA',     c('Mkt.RF', 'CMA', 'SMB', 'Mom'))
do_best_cdb(MODEL_NAME, '4F_CMA_RMW', c('Mkt.RF', 'CMA', 'SMB', 'Mom', 'RMW'))

# EW CDB -----------------------------------------------------------------

MODEL_NAME <- 'full_dynamic_std_10000'

do_ew_cdb(MODEL_NAME, '5F',          c('Mkt.RF', 'HML', 'SMB', 'RMW', 'CMA'))
do_ew_cdb(MODEL_NAME, '5F_EXCL_HML', c('Mkt.RF',        'SMB', 'RMW', 'CMA'))
do_ew_cdb(MODEL_NAME, '5F_EXCL_CMA', c('Mkt.RF', 'HML', 'SMB', 'RMW'       ))
do_ew_cdb(MODEL_NAME, '5F_EXCL_RMW', c('Mkt.RF', 'HML', 'SMB',        'CMA'))

do_ew_cdb(MODEL_NAME, '6F',          c('Mkt.RF', 'HML', 'SMB', 'Mom', 'RMW', 'CMA'))
do_ew_cdb(MODEL_NAME, '6F_EXCL_HML', c('Mkt.RF',        'SMB', 'Mom', 'RMW', 'CMA'))
do_ew_cdb(MODEL_NAME, '6F_EXCL_CMA', c('Mkt.RF', 'HML', 'SMB', 'Mom', 'RMW'       ))
do_ew_cdb(MODEL_NAME, '6F_EXCL_RMW', c('Mkt.RF', 'HML', 'SMB', 'Mom',        'CMA'))

do_ew_cdb(MODEL_NAME, '3F_HML',     c('Mkt.RF', 'HML', 'SMB'))
do_ew_cdb(MODEL_NAME, '3F_HML_RMA', c('Mkt.RF', 'HML', 'SMB', 'RMW'))
do_ew_cdb(MODEL_NAME, '3F_CMA',     c('Mkt.RF', 'CMA', 'SMB'))
do_ew_cdb(MODEL_NAME, '3F_CMA_RMW', c('Mkt.RF', 'CMA', 'SMB', 'RMW'))

do_ew_cdb(MODEL_NAME, '4F_HML',     c('Mkt.RF', 'HML', 'SMB', 'Mom'))
do_ew_cdb(MODEL_NAME, '4F_HML_RMW', c('Mkt.RF', 'HML', 'SMB', 'Mom', 'RMW'))
do_ew_cdb(MODEL_NAME, '4F_CMA',     c('Mkt.RF', 'CMA', 'SMB', 'Mom'))
do_ew_cdb(MODEL_NAME, '4F_CMA_RMW', c('Mkt.RF', 'CMA', 'SMB', 'Mom', 'RMW'))
