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

MODEL_NAME = 'oos_constant_norm_100000'

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
#' @param df.realized data frame like df.estim, including dates and 6 factors
#' of log returns, at least as long as the out-of-sample period which is given by
#' the length of the distribution or mu and sigma
#' 
#' @return saves a list of results to RData files mv_results_...
#' 
do_optimize_mv <- function(model_name = NULL, mu = NULL, sigma = NULL, strategy, selectors, df.realized) {
  # Check for correct input
  if(is.null(model_name) && any(is.null(mu), is.null(sigma))) {
    stop('No copula model passed - then optimization requires both mu and sigma')
  }
  # Whether it's on model or assumptions
  if(!is.null(model_name)) {
    based_on <- model_name
    # Load distribution simulated data and change to simple returns
    load(sprintf('data/derived/distributions/%s.RData', model_name))
    distribution_simple <- exp(distribution) - 1
    rm(distribution)
    # Subset the data using selectors
    colnames(distribution_simple) <- c('Mkt.RF', 'HML', 'SMB', 'Mom', 'RMW', 'CMA')
    distribution_simple <- distribution_simple[, selectors, ]
    # Get MV optimal weights and SRs
    mv_results <- optimize_mv(distribution = distribution_simple)
    rm(distribution_simple)
  } else {
    based_on <- 'assumption'
    # Use mu and sigma to get the mv_results
    mv_results <- optimize_mv(mu = mu, sigma = sigma)
  }
  
  colnames(mv_results$weights) <- selectors
  
  T = nrow(mv_results$weights)
  N = ncol(mv_results$weights)
  
  # Use realized returns and dates
  dates <- tail(df.realized[,'Date'], T)
  # Subset realized returns using selectors
  realized <- tail(df.realized[, selectors], T)
  # Change daily returns to simple returns
  realized <- exp(realized) - 1
  
  # Calculate the portfolios realized return
  portfolio_return <- rowSums(mv_results$weights * realized)
  
  # Save list of data points as .RData file
  results <- list(
    Date = dates$Date,
    sr = mv_results$sr,
    weights = mv_results$weights,
    portfolio_return = portfolio_return,
    based_on = based_on,
    strategy = strategy
  )
  # Give this a specific name
  save(results, file = sprintf('data/derived/mv/results_%s_%s.Rdata', based_on, strategy))
  return(results)
}

#' Takes three-dimensional objects to produce the optimal weights over time
#' Either input a distribution array 3d, or mu AND sigma in 3d
#' Returns optimized sharpe ratios and weights
#' 
#' @param distribution dim1 = Simruns, dim2 = N assets, dim3 = T periods
#' @param mu TxN vector of expected returns
#' @param sigma NxNxT array of variance-covariance matrices
#' 
#' @return list of weights and SRs
optimize_mv <- function(distribution = NULL, mu = NULL, sigma = NULL) {
  # Check for correct input
  if(is.null(distribution) && any(is.null(mu), is.null(sigma))) {
    stop('No distribution passed - then optimization requires both mu and sigma')
  }
  if(!is.null(distribution)) {
    # If distribution is passed
    times <- 1:dim(distribution)[3]
    N <- ncol(distribution)
  } else {
    # If mu and sigma are passed
    times <- 1:nrow(mu)
    N <- ncol(mu)
  }
  
  T <- length(times)
  weights <- matrix(NA, ncol = N, nrow = T)
  sr <- rep(NA, T)
  
  if(!is.null(distribution)) {
    # If distribution is passed
    for (t in times) {
      tic(sprintf('Optimal weights at t = %d', t))
      op <- optimize_mv_1_period(distribution_t = distribution[,,t])
      toc()
      
      if (op$convergence > 0) {
        warning(sprintf('No convergence for t = %d', t))
      }
      
      weights[t, ] <- op$pars
      sr[t] <- -tail(op$values,1)
    }
  } else {
    # If mu and sigma are passed
    for (t in times) {
      tic(sprintf('Optimal weights at t = %d', t))
      op <- optimize_mv_1_period(mu_t = mu[t,], sigma_t = sigma[,,t])
      toc()
      
      if (op$convergence > 0) {
        warning(sprintf('No convergence for t = %d', t))
      }
      
      weights[t, ] <- op$pars
      sr[t] <- -tail(op$values,1)
    }
  }
  
  # Collect in list
  list(
    weights = weights, sr = sr
  )
}

#' Returns optimized weights from solnp, maximizing SR
#' Can take either a matrix of simulated data or
#' a given sigma and mu.
#' Is constrained, weights sum to 1 and all greater or equal to zero
#' 
#' @param distribution_t (SimRuns x N) matrix of simulated data (1 period)
#' @param mu_t N length vector of expected returns (simple)
#' @param sigma_t NxN variance-covariance matrix
#'
#' @return weights optimized for highest SR
optimize_mv_1_period <- function(distribution_t = NULL, mu_t = NULL, sigma_t = NULL, eqfun = sum, eqB = 1,
                                 LB = NULL, x0 = NULL, ...) {
  # Check if simulated data inputted or sigma, mu
  if (is.null(distribution_t)) {
    # Mu and sigma as inputted
    N = length(mu_t)
  } else {
    # Inputs needed from simulated 1-step-ahead distribution
    mu_t <- colMeans(distribution_t)
    sigma_t <- cov(distribution_t)
    N <- ncol(distribution_t)  
  }
  
  # No negative weights
  if (is.null(LB)) {
    LB <- rep(0, N)
  }
  
  # Start with EW portfolio
  if (is.null(x0)) {
    x0 <- rep(1 / N, N)
  }
  
  fn <- function(weights, mu_t, sigma_t) -sharpe_ratio(weights, mu_t, sigma_t)[1]
  
  solnp(
    x0,
    fn,
    eqfun = eqfun,
    eqB = eqB,
    LB = LB,
    
    # Optimizer keeps climbing the multidimensional mountains with tiny
    # tiny steps. Take some strides (default delta = 1e-7)!!!
    control = list(trace = 0,
                   delta = 1e-5),
    mu_t = mu_t,
    sigma_t = sigma_t
  )
}

#' Calculates sharpe ratio
#' 
#' @param weights N length vector of weights
#' @param mu_t N length vector of expected returns (simple)
#' @param sigma_t NxN variance-covariance matrix
#' 
#' @return sr scalar - sharpe ratio
sharpe_ratio <- function(weights, mu_t, sigma_t) {
  # 
  sr <- (weights %*% mu_t) / sqrt(weights %*% sigma_t %*% weights)
}

#' Get the fixed weights portfolio results for comparison
#' with MV optimized portfolios. E.g. equal weights..
#' 
#' Returns a list with dates, chosen fixed weights, and 
#' resulting realized portfolio return (based on df.realized's
#' last T values)
#' 
#' @param model_name copula model, one of 'ghskt', 'ght', 'norm'
#' @param strategy name for out file for this strategy
#' @param selectors chosen allowed assets in this strategy, 
#' choosing from Mkt.RF, HML, SMB, Mom, RMW, CMA
#' @param weights vector of length 6 (all assets) with weights, sum to 1
#' and all greater or equal to zero
#' @param T length of out-of-sample period
#' @param df.realized data frame like df.estim, including dates and 6 factors
#' of log returns, at least as long as the out-of-sample period T
#' 
#' @return saves a list of results to RData files mv_results_...
#' 
do_fixed_weights_mv <- function(strategy, selectors, weights_fixed, T, df.realized) {
  # Error check on weights
  if(sum(weights_fixed) != 1 | any(weights_fixed < 0)) {
    stop('Weights incorrectly inputted - need to sum to one and all be greater or equal to zero')
  }
  # Preliminary constants. Count number of active assets in selectors
  based_on <- 'fixed_weights'
  N_total = ncol(df.realized[,-1])
  N_active = length(selectors)
  
  # Create results with weights. Sharpe ratio from optimization undefined here
  mv_results <- list(
    weights = matrix(weights_fixed, ncol = N_total, nrow = T),
    sr = rep(NA, T)
  )
  # Name assets
  colnames(mv_results$weights) <- colnames(df.realized[,-1])
  
  # Realized returns and dates
  dates <- tail(df.realized[,'Date'], T)
  # Subset realized returns using selectors
  realized <- tail(df.realized[, selectors], T)
  # Change daily returns to simple returns
  realized <- exp(realized) - 1
  
  # Calculate the portfolios realized return
  portfolio_return <- rowSums(mv_results$weights * realized)
  
  # Save list of data points as .RData file
  results <- list(
    Date = dates$Date,
    sr = mv_results$sr,
    weights = mv_results$weights,
    portfolio_return = portfolio_return,
    based_on = based_on,
    strategy = strategy
  )
  # Give this a specific name
  save(results, file = sprintf('data/derived/mv/results_%s_%s.Rdata', based_on, strategy))
  return(results)
}

# Get optimization results for different copula portfolios ------------------------------------

load('data/derived/weekly-estim.RData')
results_All <- do_optimize_mv(model_name = MODEL_NAME, strategy = 'All', selectors = c('Mkt.RF','HML','SMB','Mom','RMW','CMA'), df.realized = df.estim)
results_Four_HML <- do_optimize_mv(model_name = MODEL_NAME, strategy = 'Four+HML', selectors = c('Mkt.RF', 'HML', 'SMB', 'Mom', 'RMW'), df.realized = df.estim)
results_Four_CMA <- do_optimize_mv(model_name = MODEL_NAME, strategy = 'Four+CMA', selectors = c('Mkt.RF', 'SMB', 'Mom', 'RMW', 'CMA'), df.realized = df.estim)

results_5F_HML <- do_optimize_mv(model_name = MODEL_NAME, strategy = 'FF_5F_INCL_HML', selectors = c('Mkt.RF','HML','SMB','RMW','CMA'), df.realized = df.estim)
results_5F <- do_optimize_mv(model_name = MODEL_NAME, strategy = 'FF_5F_EXCL_HML', selectors = c('Mkt.RF','SMB','RMW','CMA'), df.realized = df.estim)

# And for EW portfolio ----

results_EW_6F <- do_fixed_weights_mv(strategy = 'EW_6F', T = 2765, 
                                  selectors = c('Mkt.RF','HML','SMB','Mom','RMW','CMA'), 
                                  weights_fixed = rep(1/6, 6),
                                  df.realized = df.estim
                                  )

# Get optimization results for different assumed mu and sigma --------------------------
load('data/derived/weekly-estim.RData')

# Function to create mu and sigma from sample

.sample_mu <- function(logreturndata, nOOS) {
  ret <- exp(logreturndata) - 1
  mu <- colMeans(ret)
  matrix(mu, nrow = nOOS, ncol = length(mu), byrow = TRUE)
}

.sample_sigma <- function(logreturndata, nOOS) {
  ret <- exp(logreturndata) - 1
  sigma = cov(ret)
  array(sigma, dim = c(ncol(sigma), ncol(sigma), nOOS))
}

results_Sample <- do_optimize_mv(mu = .sample_mu(df.estim[,-1], 2765), sigma = .sample_sigma(df.estim[,-1], 2765), strategy = 'Sample',
               selectors = c('Mkt.RF','HML','SMB','Mom','RMW','CMA'), df.realized = df.estim)
results_Sample_5F <- do_optimize_mv(mu = .sample_mu(df.estim[,c(-1,-5)], 2765), sigma = .sample_sigma(df.estim[,c(-1,-5)], 2765), strategy = 'Sample',
                                     selectors = c('Mkt.RF','HML','SMB','RMW','CMA'), df.realized = df.estim)


# Load and use FF MV results ----------------------------------------------
load('data/derived/mv_plot_data.RData')

load_results <- function(name) {
  load(sprintf('data/derived/mv/results_oos_%s.RData', name))
  results
}

result_names <- list(
  'indep_100000_All',
  'constant_norm_10000_All',
  'constant_std_10000_All',
  'constant_ghst_10000_All',
  'dynamic_norm_10000_All',
  'dynamic_std_10000_All'
)

results <- lapply(names, load_results)
names(results) <- result_names

# Weights over time, smoothed and unsmoothed ------------------------------

.weights_smooth_plot <- function(results1, results2, name1, name2) {
  
  prepare_plotdf <- function(results, name) {
    dates <- results$Date
    results <- gather(as.data.frame(results$weights), 'factor','weight')
    results$Date <- rep(dates, length(unique(results$factor)))
    results$factor <- factor(results$factor, levels = c('Mkt.RF','HML','SMB','Mom','RMW','CMA'))
    results$id <- name
    return(results)
  }
  
  results1 <- prepare_plotdf(results1, name1)
  results2 <- prepare_plotdf(results2, name2)

  g <- ggplot(mapping = aes(x = Date, y = weight, color = id))+
    geom_smooth(data = results1)+
    geom_smooth(data = results2)+
    theme_Publication()+
    scale_colour_Publication()+
    coord_cartesian(ylim = c(0, 1))+
    facet_wrap(~ factor, nrow = 2, ncol = 3)
  
  
}

.weights_plot <- function(results1, results2, name1, name2) {
  
  prepare_plotdf <- function(results, name) {
    dates <- results$Date
    results <- gather(as.data.frame(results$weights), 'factor','weight')
    results$Date <- rep(dates, length(unique(results$factor)))
    results$factor <- factor(results$factor, levels = c('Mkt.RF','HML','SMB','Mom','RMW','CMA'))
    results$id <- name
    return(results)
  }
  
  results1 <- prepare_plotdf(results1, name1)
  results2 <- prepare_plotdf(results2, name2)
  
  g <- ggplot(mapping = aes(x = Date, y = weight, color = id))+
    geom_line(data = results1)+
    geom_line(data = results2)+
    theme_Publication()+
    scale_colour_Publication()+
    coord_cartesian(ylim = c(0, 1))+
    facet_wrap(~ factor, nrow = 2, ncol = 3)
  
  
}
g <- .weights_plot(results_5F_HML, results_Sample_5F, 'In-sample simulated dynamic ghskt copula','In-sample sample mu and sigma')
g <- .weights_smooth_plot(results_5F_HML, results_Sample_5F, 'In-sample simulated dynamic ghskt copula','In-sample sample mu and sigma')

g <- .weights_plot(results_Four_HML,results_Four_CMA, 'HML', 'CMA')
g <- .weights_smooth_plot(results_Four_HML,results_Four_CMA, 'Five factors HML', 'Five factors CMA')

g <- .weights_plot(results_Sample, results_All, 'Sample', 'Simulated')
g <- .weights_smooth_plot(results_Sample, results_All, 'Sample', 'Simulated')


# Return over time --------------------------------------------------------
#Something is funky here
.cumret_plot <- function(results1, results2, name1, name2) {
  
  prepare_plotdf <- function(results, name) {
    results <- data.frame(ret = results$portfolio_return,
                          cumret = c(1,
                                     cumprod(
                                       1+results$portfolio_return[2:length(results$portfolio_return)]
                                       )
                                     ),
                          id = name,
                          Date = results$Date)
    return(results)
  }
  
  results1 <- prepare_plotdf(results1, name1)
  results2 <- prepare_plotdf(results2, name2)
  
  g <- ggplot(mapping = aes(x = Date, y = cumret, color = id))+
    geom_smooth(data = results1)+
    geom_smooth(data = results2)+
    theme_Publication()+
    scale_colour_Publication()
  
  return(g)
}

g <- .cumret_plot(results_All, results_Sample, 'Simulated', 'Sample')
g <- .cumret_plot(results_Sample, results_EW_6F, 'Sample', 'Equal weighted')
#  Density plot ------------------------------------------------------------------------
.density_plot <- function(results1, results2, name1, name2) {
  
  prepare_plotdf <- function(results, name) {
    results <- data.frame(ret = results$portfolio_return,
                          id = name,
                          Date = results$Date)
    return(results)
  }
  
  results1 <- prepare_plotdf(results1, name1)
  results2 <- prepare_plotdf(results2, name2)
  
  g <- ggplot(mapping = aes(ret, colour = id))+
    geom_density(alpha = 0.1, data = results1)+
    geom_density(alpha = 0.1, data = results2)+
    theme_Publication()+
    scale_colour_Publication()
    
}

g <- .density_plot(results_5F, results_5F_HML, '5F', '5F+HML')
g <- .density_plot(results_Sample, results_All, 'Sample', 'Simulated')
g <- .density_plot(results_Four_CMA, results_Four_HML, 'Four CMA', 'Four HML')
g <- .density_plot(results_Sample, results_Sample_5F, 'Sample All', 'Sample ex Momentum')
g <- .density_plot(results_Sample, results_EW_6F, 'Sample All','Equal weighted All')
#  ------------------------------------------------------------------------
# Table of summary stats, MDD of returns, realized SR, mean return, sd, mean weights etc

.summary_stats <- function(result) {
  ret <- result$portfolio_return
  # get the standard error on skewness, kurtosis?
  stats_list <- c('nobs','Maximum','Minimum','Mean','Median','Stdev','Skewness','Kurtosis')
  
  table <- ret %>%
    basicStats() %>%
      .[stats_list,] %>%
        round(., digits = 4)
  names(table) <- stats_list
  # add maximum drawdown
  # box plot, percentile ranges, oos
  out <- as.data.frame(t(c(table, SR = (52*mean(ret)) / (sqrt(52)*sd(ret)), colMeans(result$weights))))
  
}

# results_list = list(
#   results_EW_6F = results_EW_6F,
#   results_All = results_All,
#   results_Sample = results_Sample,
#   results_5F_HML = results_5F_HML,
#   results_Sample_5F = results_Sample_5F,
#   results_Four_CMA = results_Four_CMA,
#   results_Four_HML = results_Four_HML,
#   results_5F = results_5F
# )

summary_table <- bind_rows(lapply(results, .summary_stats), .id = 'id')
stargazer(summary_table, type = 'text', summary = FALSE)
