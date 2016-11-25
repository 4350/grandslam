#' What we need to do here
#' Simulate for all copulae, and maybe for independence copula, as well as garchs by themselves.
#' Do this simulation out of sample, leaving e.g. everything after 2000 as the OOS
#' 
#' Cumulative returns of simulated stuff looks crazy high - what's going on?
#' How would we impose restrictions on how fast weights can change?
#' moving average? can you rebalance any amount? yes, the amount should be about depth in the market,
#' can you sell the whole fund without price impact, not about number of trades, will be the same
#' even if we only rebalance 1%. so.. the rebalancing limits is rather about getting robust
#' performance, following equal weights mindset, risk parity etc. then maybe 2x equal weights? 50%?
#' Think about how to scale sharpe ratio. Should be oK as we do.

# Library and setup -------------------------------------------------------
rm(list = ls())
library(tictoc)
library(Rsolnp)
library(fBasics)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggfortify)
library(gridExtra)
library(devtools)
library(extrafont)
library(PerformanceAnalytics)
library(stargazer)

load_all('wimbledon')
load_all('australian')

load('data/derived/weekly-estim.RData')
MODEL_NAME = 'full_dynamic_std_10000'

# Functions used in optimization  ---------------------------------------------------------------

#' Get the MV optimal portfolio results for a copula simulated 
#' distribution, or for assumptions of mu and sigma over time.
#' 
#' Returns a list with dates, sharpe ratios (maximized by the optimizer),
#' optimal weights constrained to sum to 1 and be greater or equal
#' to zero, and the resulting realized portfolio return (based on df.realized's
#' last T values)
#' 
#' @param model_name copula model, one of 'ghskt', 'ght', 'norm'
#' @param strategy name for out file for this strategy
#' @param selectors chosen allowed assets in this strategy, 
#' choosing from Mkt.RF, HML, SMB, Mom, RMW, CMA
#' @param q CDB VaR cutoff
#' 
#' @return saves a list of results to RData files mv_results_...
#' 
do_optimize_mv <- function(model_name, strategy, selectors, q = 0.05) {
  
  # Load distribution simulated data
  load(sprintf('data/derived/distributions/%s.RData', model_name))
  
  # Subset the data using selectors
  colnames(distribution) <- c('Mkt.RF', 'HML', 'SMB', 'Mom', 'RMW', 'CMA')
  distribution <- distribution[, selectors, ]
  
  times <- 1:dim(distribution)[3]
  N <- ncol(distribution)
  
  T <- length(times)
  
  # Outputs
  weights <- matrix(NA, ncol = N, nrow = T)
  sr <- rep(NA, T)
  cdb <- rep(NA, T)
  
  # Get MV optimal weights and SRs and CDBs
  
  for (t in times) {
    tic(sprintf('Optimal weights at t = %d', t))
    
    # Mu, sigma from dist
    mu_t <- colMeans(distribution[,,t])
    sigma_t <- cov(distribution[,,t])
    
    # Run optimizer
    op <- optimize_mv_constrOptim(weights = rep(1/N,N), fn = sharpe_ratio_fn(mu_t, sigma_t))
    
    toc()
    
    # Get CDB
    
    fn_cdb <- cdb_fn(q, distribution[,,t])
    cdb[t] <- fn_cdb(op$weights)
    
    # Save this
    
    weights[t, ] <- op$weights
    sr[t] <- op$sr
  }
  
  colnames(weights) <- selectors
  
  # Get realized dates, returns and CDB
  metrics <- portfolio_metrics(weights, distribution, q, selectors)
  
  # Save list of data points as .RData file
  results <- c(
    metrics,
    list(
      sr_optim = sr,
      model_name = model_name,
      strategy = strategy,
      weights = weights
    )
  )
  
  # Give this a specific name
  save(results, file = sprintf('data/derived/mv/results_%s_%s.Rdata', model_name, strategy))
  return(results)
}

optimize_mv_constrOptim <- function(weights, fn) {
  N <- length(weights)
  
  # Don't optimize on the first weight, as it is restricted by the
  # sum(weights) == 1 condition.
  weights1 <- weights[-1]
  
  # Call constrOptim with different methods
  optimize <- function(weights1, method = 'Nelder-Mead') {
    constrOptim(
      weights1,
      function(w) fn(c(1 - sum(w), w)),
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
    sr = -op$value
  )
}

#' Calculates sharpe ratio
#' 
#' @param mu_t N length vector of expected excess returns
#' @param sigma_t NxN variance-covariance matrix
#' 
#' @return function for optimization, negative sharpe ratio
sharpe_ratio_fn <- function(mu_t, sigma_t) {
  # 
  function(weights) {
    sr <- (weights %*% mu_t) / sqrt(weights %*% sigma_t %*% weights)
    -sr
  }
  
}

#' Optimize with fixed mu, sigma
#' 
#' @param model_name copula model, one of 'ghskt', 'ght', 'norm'
#' @param strategy name for out file for this strategy
#' @param selectors chosen allowed assets in this strategy, 
#' choosing from Mkt.RF, HML, SMB, Mom, RMW, CMA
#' @param df.sample sample data frame (df.estim)
#' @param q CDB VaR cutoff

do_optimize_fixed <- function(model_name, strategy, selectors, df.sample = df.estim, q = 0.05) {
  
  # Load distribution simulated data
  load(sprintf('data/derived/distributions/%s.RData', model_name))
  
  # Subset the data using selectors
  colnames(distribution) <- c('Mkt.RF', 'HML', 'SMB', 'Mom', 'RMW', 'CMA')
  distribution <- distribution[, selectors, ]
  
  # Select subset and shorten to tail_T length
  df.sample <- df.sample[ , selectors]
  
  N <- ncol(df.sample)
  T <- nrow(df.sample)-1 #adj due to one day falls off in distributions
  
  # Get MV optimal weights and SRs
  mu = colMeans(df.sample)
  sigma = cov(df.sample)
  op <- optimize_mv_constrOptim(weights = rep(1/N,N), fn = sharpe_ratio_fn(mu, sigma))
    
  # Outputs
  weights <- matrix(op$weights, ncol = N, nrow = T, byrow = TRUE)
  sr <- rep(op$sr, T)
  
  colnames(weights) <- selectors
  
  # Get dates realized returns and CDB
  metrics <- portfolio_metrics(weights, distribution, q, selectors)
  
  # Save list of data points as .RData file
  results <- c(
    metrics,
    list(
      sr_optim = sr,
      model_name = model_name,
      strategy = strategy,
      weights = weights
    )
  )
  
  # Give this a specific name
  save(results, file = sprintf('data/derived/mv/results_%s_%s.Rdata', 'sample', strategy))
  return(results)
}

# Get optimization results dynamic std full copula means ------------------------------------

# Without Momentum
do_optimize_mv(MODEL_NAME, '5F'         , c('Mkt.RF', 'HML', 'SMB', 'RMW', 'CMA'))
do_optimize_mv(MODEL_NAME, '5F_EXCL_HML', c('Mkt.RF',        'SMB', 'RMW', 'CMA'))
do_optimize_mv(MODEL_NAME, '5F_EXCL_CMA', c('Mkt.RF', 'HML', 'SMB', 'RMW'       ))
do_optimize_mv(MODEL_NAME, '5F_EXCL_RMW', c('Mkt.RF', 'HML', 'SMB',        'CMA'))

# With Momentum
do_optimize_mv(MODEL_NAME, '6F'         , c('Mkt.RF', 'HML', 'SMB', 'Mom', 'RMW', 'CMA'))
do_optimize_mv(MODEL_NAME, '6F_EXCL_HML', c('Mkt.RF',        'SMB', 'Mom', 'RMW', 'CMA'))
do_optimize_mv(MODEL_NAME, '6F_EXCL_CMA', c('Mkt.RF', 'HML', 'SMB', 'Mom', 'RMW'       ))
do_optimize_mv(MODEL_NAME, '6F_EXCL_RMW', c('Mkt.RF', 'HML', 'SMB', 'Mom',        'CMA'))

# Optimize sample full 5F ---------------------------------------------------------

do_optimize_fixed(MODEL_NAME, '5F'         , c('Mkt.RF', 'HML', 'SMB', 'RMW', 'CMA'))
do_optimize_fixed(MODEL_NAME, '5F_EXCL_HML', c('Mkt.RF',        'SMB', 'RMW', 'CMA'))
do_optimize_fixed(MODEL_NAME, '5F_EXCL_CMA', c('Mkt.RF', 'HML', 'SMB', 'RMW'       ))
do_optimize_fixed(MODEL_NAME, '5F_EXCL_RMW', c('Mkt.RF', 'HML', 'SMB',        'CMA'))


# Optimize sample full 6F ---------------------------------------------------------

do_optimize_fixed(MODEL_NAME, '6F'         , c('Mkt.RF', 'HML', 'SMB', 'Mom', 'RMW', 'CMA'))
do_optimize_fixed(MODEL_NAME, '6F_EXCL_HML', c('Mkt.RF',        'SMB', 'Mom', 'RMW', 'CMA'))
do_optimize_fixed(MODEL_NAME, '6F_EXCL_CMA', c('Mkt.RF', 'HML', 'SMB', 'Mom', 'RMW'       ))
do_optimize_fixed(MODEL_NAME, '6F_EXCL_RMW', c('Mkt.RF', 'HML', 'SMB', 'Mom',        'CMA'))



