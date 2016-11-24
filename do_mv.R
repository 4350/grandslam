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
#' @param mu TxN vector of expected returns
#' @param sigma NxNxT array of variance-covariance matrices
#' @param strategy name for out file for this strategy
#' @param selectors chosen allowed assets in this strategy, 
#' choosing from Mkt.RF, HML, SMB, Mom, RMW, CMA
#' @param df.realized data frame like df.estim, including dates and 6 factors
#' of log returns, at least as long as the out-of-sample period which is given by
#' the length of the distribution or mu and sigma
#' @param df.mean optional df of factors used to calculate mean if 
#' copula estimate from mean equation is not to be used
#' 
#' @return saves a list of results to RData files mv_results_...
#' 
do_optimize_mv <- function(model_name, strategy, selectors, df.realized, df.mean = NULL) {
  
  # Load distribution simulated data
  load(sprintf('data/derived/distributions/%s.RData', model_name))
  
  # Simplify df.mean if inputted. Select selectors
  if(!is.null(df.mean)) {
    df.mean <- df.mean[,-1]
    df.mean <- df.mean[, selectors]
  }
  
  # Subset the data using selectors
  colnames(distribution) <- c('Mkt.RF', 'HML', 'SMB', 'Mom', 'RMW', 'CMA')
  distribution <- distribution[, selectors, ]
  
  times <- 1:dim(distribution)[3]
  N <- ncol(distribution)
  
  T <- length(times)
  # Outputs
  weights <- matrix(NA, ncol = N, nrow = T)
  sr <- rep(NA, T)
  
  # Get MV optimal weights and SRs
  
  for (t in times) {
    tic(sprintf('Optimal weights at t = %d', t))
    
    # Mu either from mean or copula distribution
    
    if(!is.null(df.mean)) {
      mu_t <- colMeans(df.mean)
    } else {
      mu_t <- colMeans(distribution[,,t])
    }
    
    # Sigma always from copula distribution
    
    sigma_t <- cov(distribution[,,t])
    
    # Run optimizer
    op <- optimize_mv_constrOptim(weights = rep(1/N,N), fn = sharpe_ratio_fn(mu_t, sigma_t))
    
    toc()
    
    weights[t, ] <- op$weights
    sr[t] <- op$sr
  }
  
  colnames(weights) <- selectors
  
  # Use realized returns and dates
  dates <- tail(df.realized[,'Date'], T)
  # Subset realized returns using selectors
  realized <- tail(df.realized[, selectors], T)
  
  # Calculate the portfolios realized return
  portfolio_return <- rowSums(weights * realized)
  
  # Save list of data points as .RData file
  results <- list(
    Date = dates$Date,
    sr = sr,
    weights = weights,
    portfolio_return = portfolio_return,
    based_on = model_name,
    strategy = strategy
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
do_optimize_fixed <- function(model_name, strategy, selectors, df.sample, df.realized) {
  
  # Save the dates and then select
  dates <- df.realized[,'Date']
  # Select returns 
  df.realized <- df.realized[ , selectors]
  df.sample <- df.sample[ , selectors]
  
  N <- ncol(df.realized)
  T <- nrow(df.realized)
  
  # Get MV optimal weights and SRs
  mu = colMeans(df.sample)
  sigma = cov(df.sample)
  op <- optimize_mv_constrOptim(weights = rep(1/N,N), fn = sharpe_ratio_fn(mu, sigma))
    
  # Outputs
  weights <- matrix(op$weights, ncol = N, nrow = T, byrow = TRUE)
  sr <- rep(op$sr, T)
  
  colnames(weights) <- selectors
  
  # Calculate the portfolios realized return
  portfolio_return <- rowSums(weights * df.realized)
  
  # Save list of data points as .RData file
  results <- list(
    Date = dates$Date,
    sr = sr,
    weights = weights,
    portfolio_return = portfolio_return,
    based_on = model_name,
    strategy = strategy
  )
  # Give this a specific name
  save(results, file = sprintf('data/derived/mv/results_%s_%s.Rdata', model_name, strategy))
  return(results)
}

# Get optimization results dynamic std full copula means ------------------------------------

load('data/derived/weekly-estim.RData')
MODEL_NAME = 'full_dynamic_std_10000'

# Without Momentum
do_optimize_mv(MODEL_NAME, strategy = '5F',
               selectors = c('Mkt.RF', 'HML', 'SMB', 'RMW', 'CMA'),
               df.realized = df.estim)

do_optimize_mv(MODEL_NAME, strategy = '5F_EXCL_HML',
               selectors = c('Mkt.RF',        'SMB', 'RMW', 'CMA'),
               df.realized = df.estim)

do_optimize_mv(MODEL_NAME, strategy = '5F_EXCL_CMA',
               selectors = c('Mkt.RF', 'HML', 'SMB', 'RMW'),
               df.realized = df.estim)

do_optimize_mv(MODEL_NAME, strategy = '5F_EXCL_RMW',
               selectors = c('Mkt.RF', 'HML', 'SMB', 'CMA'),
               df.realized = df.estim)

# With Momentum
do_optimize_mv(MODEL_NAME, strategy = '6F',
               selectors = c('Mkt.RF', 'HML', 'SMB', 'Mom', 'RMW', 'CMA'),
               df.realized = df.estim)

do_optimize_mv(MODEL_NAME, strategy = '6F_EXCL_HML',
               selectors = c('Mkt.RF', 'SMB', 'Mom', 'RMW', 'CMA'),
               df.realized = df.estim)

do_optimize_mv(MODEL_NAME, strategy = '6F_EXCL_CMA',
               selectors = c('Mkt.RF', 'HML', 'SMB', 'Mom', 'RMW'),
               df.realized = df.estim)

do_optimize_mv(MODEL_NAME, strategy = '6F_EXCL_RMW',
               selectors = c('Mkt.RF', 'HML', 'SMB', 'Mom', 'CMA'),
               df.realized = df.estim)

# Optimize sample full 5F ---------------------------------------------------------

load('data/derived/weekly-estim.RData')
do_optimize_fixed('full_sample', strategy = '5F',
                  selectors = c('Mkt.RF', 'HML', 'SMB', 'RMW', 'CMA'),
                  df.sample = df.estim[1:2766,],
                  df.realized = df.estim[2:2766,])

do_optimize_fixed('full_sample', strategy = '5F_EXCL_HML',
                  selectors = c('Mkt.RF',        'SMB', 'RMW', 'CMA'),
                  df.sample = df.estim[1:2766,],
                  df.realized = df.estim[2:2766,])

do_optimize_fixed('full_sample', strategy = '5F_EXCL_CMA',
                  selectors = c('Mkt.RF', 'HML', 'SMB', 'RMW'),
                  df.sample = df.estim[1:2766,],
                  df.realized = df.estim[2:2766,])

do_optimize_fixed('full_sample', strategy = '5F_EXCL_RMW',
                  selectors = c('Mkt.RF', 'HML', 'SMB', 'CMA'),
                  df.sample = df.estim[1:2766,],
                  df.realized = df.estim[2:2766,])



# Optimize sample full 6F ---------------------------------------------------------

load('data/derived/weekly-estim.RData')
do_optimize_fixed('full_sample', strategy = '6F',
                  selectors = c('Mkt.RF', 'HML', 'SMB', 'Mom', 'RMW', 'CMA'),
                  df.sample = df.estim[1:2766,],
                  df.realized = df.estim[2:2766,])

do_optimize_fixed('full_sample', strategy = '6F_EXCL_HML',
                  selectors = c('Mkt.RF', 'SMB', 'Mom', 'RMW', 'CMA'),
                  df.sample = df.estim[1:2766,],
                  df.realized = df.estim[2:2766,])

do_optimize_fixed('full_sample', strategy = '6F_EXCL_CMA',
                  selectors = c('Mkt.RF', 'HML', 'SMB', 'Mom', 'RMW'),
                  df.sample = df.estim[1:2766,],
                  df.realized = df.estim[2:2766,])

do_optimize_fixed('full_sample', strategy = '6F_EXCL_RMW',
                  selectors = c('Mkt.RF', 'HML', 'SMB', 'Mom', 'CMA'),
                  df.sample = df.estim[1:2766,],
                  df.realized = df.estim[2:2766,])



