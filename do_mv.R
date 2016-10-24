  
  
# Library and setup -------------------------------------------------------
rm(list = ls())
library(tictoc)
library(Rsolnp)

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggfortify)
library(gridExtra)
library(devtools)
library(extrafont)

load_all('wimbledon')

MODEL_NAME = 'dynamic_ghskt'

# Functions ---------------------------------------------------------------

#' Get the MV optimal portfolio results for a copula simulated 
#' distribution, or for assumptions of mu and sigma over time.
#' 
#' Returns a list with dates, sharpe ratios (maximized by the optimizer),
#' optimal weights constrained to sum to 1 and be greater or equal
#' to zero, and the resulting realized portfolio return (based on df.estim)
#' 
#' @param model_name copula model, one of 'ghskt', 'ght', 'norm'
#' @param strategy name for out file for this strategy
#' @param selectors chosen allowed assets in this strategy, 
#' choosing from Mkt.RF, HML, SMB, Mom, RMW, CMA
#' 
#' @return saves a list of results to RData files mv_results_...
#' 
do_optimize_mv <- function(model_name = NULL, mu = NULL, sigma = NULL, strategy, selectors) {
  # Check for correct input
  if(is.null(model_name) && any(is.null(mu), is.null(sigma))) {
    stop('No copula model passed - then optimization requires both mu and sigma')
  }
  if(!is.null(model_name)) {
    based_on <- model_name
    # Load distribution simulated data and change to simple returns
    load(sprintf('data/derived/distribution_%s.RData', model_name))
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
  
  # Load realized returns and dates
  load('data/derived/weekly-estim.RData')
  dates <- tail(df.estim[,'Date'], T)
  # Subset realized returns using selectors
  realized <- tail(df.estim[, selectors], T)
  # Change daily returns to simple returns
  realized <- exp(realized) - 1
  rm(df.estim)
  
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
  save(results, file = sprintf('data/derived/mv_results_%s_%s.Rdata', based_on, strategy))
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

# Get results for different copula portfolios ------------------------------------

results_All <- do_optimize_mv(model_name = MODEL_NAME, strategy = 'All', selectors = c('Mkt.RF','HML','SMB','Mom','RMW','CMA'))
results_Four_HML <- do_optimize_mv(model_name = MODEL_NAME, strategy = 'Four+HML', selectors = c('Mkt.RF', 'HML', 'SMB', 'Mom', 'RMW'))
results_Four_CMA <- do_optimize_mv(model_name = MODEL_NAME, strategy = 'Four+CMA', selectors = c('Mkt.RF', 'SMB', 'Mom', 'RMW', 'CMA'))

results_5F_HML <- do_optimize_mv(model_name = MODEL_NAME, strategy = 'FF_5F_INCL_HML', selectors = c('Mkt.RF','HML','SMB','RMW','CMA'))
results_5F <- do_optimize_mv(model_name = MODEL_NAME, strategy = 'FF_5F_EXCL_HML', selectors = c('Mkt.RF','SMB','RMW','CMA'))


# Get results for different assumed mu and sigma --------------------------
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

results_Sample_All <- do_optimize_mv(mu = .sample_mu(df.estim[,-1], 2765), sigma = .sample_sigma(df.estim[,-1], 2765), strategy = 'Sample',
               selectors = c('Mkt.RF','HML','SMB','Mom','RMW','CMA'))
results_Sample_5F <- do_optimize_mv(mu = .sample_mu(df.estim[,c(-1,-5)], 2765), sigma = .sample_sigma(df.estim[,c(-1,-5)], 2765), strategy = 'Sample',
                                     selectors = c('Mkt.RF','HML','SMB','RMW','CMA'))

# Load and use FF MV results ----------------------------------------------

weights_smooth_plot <- function(results1, results2, name1, name2) {
  
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

weights_plot <- function(results1, results2, name1, name2) {
  
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
g <- weights_plot(results_5F_HML, results_Sample_5F, 'In-sample simulated dynamic ghskt copula','In-sample sample mu and sigma')
g <- weights_smooth_plot(results_5F_HML, results_Sample_5F, 'In-sample simulated dynamic ghskt copula','In-sample sample mu and sigma')

g <- weights_plot(results_Four_HML,results_Four_CMA, 'HML', 'CMA')
g <- weights_smooth_plot(results_Four_HML,results_Four_CMA, 'Five factors HML', 'Five factors CMA')

g <- weights_plot(results_Sample_All, results_All, 'Sample', 'Simulated')
g <- weights_smooth_plot(results_Sample_All, results_All, 'Sample', 'Simulated')
# Graph of return over time
# Graph of diff in return over time
# Graph of density of return comparison
# Table of summary stats, MDD of returns, realized SR, mean return, sd, mean weights etc


#  ------------------------------------------------------------------------



rm(list = ls())

load('data/derived/mv_results_FF_5F_EXCL_HML_dynamic_ghskt.Rdata')
FF_5F_EXCL_HML <- out.list
load('data/derived/mv_results_FF_5F_INCL_HML_dynamic_ghskt.Rdata')
FF_5F_INCL_HML <- out.list
rm(out.list)

plot(FF_5F_EXCL_HML$portfolio_return, type = 'l', col = 'red')
lines(FF_5F_INCL_HML$portfolio_return, type = 'l', col = 'blue')

plot(FF_5F_INCL_HML$portfolio_return - FF_5F_EXCL_HML$portfolio_return, type = 'l')
mean(FF_5F_INCL_HML$portfolio_return - FF_5F_EXCL_HML$portfolio_return)

d1 <- density(FF_5F_EXCL_HML$portfolio_return)
d2 <- density(FF_5F_INCL_HML$portfolio_return)
plot(range(d1$x, d2$x), range(d1$y, d2$y), type = "n", xlab = "x",
     ylab = "Density")
lines(d1, col = "red")
lines(d2, col = "blue")

plot(FF_5F_EXCL_HML$weights[,'CMA'], type = 'l', col = 'red')
lines(FF_5F_INCL_HML$weights[,'CMA'], type = 'l', col = 'blue')

plot(cumprod(1+FF_5F_EXCL_HML$portfolio_return))
plot(cumprod(1+FF_5F_INCL_HML$portfolio_return))

mean(FF_5F_EXCL_HML$portfolio_return)/sd(FF_5F_EXCL_HML$portfolio_return)
mean(FF_5F_INCL_HML$portfolio_return)/sd(FF_5F_INCL_HML$portfolio_return)


# Load and use MV results -------------------------------------------------
rm(list = ls())

load('data/derived/mv_results_All_dynamic_ghskt.Rdata')
six_all <- out.list
load('data/derived/mv_results_Four+CMA_dynamic_ghskt.Rdata')
five_CMA <- out.list
load('data/derived/mv_results_Four+HML_dynamic_ghskt.Rdata')
five_HML <- out.list

rm(out.list)

plot(five_HML$portfolio_return, type = 'l', col = 'red')
lines(five_CMA$portfolio_return, type = 'l', col = 'blue')

plot(five_HML$portfolio_return - five_CMA$portfolio_return, type = 'l')

d1 <- density(five_HML$portfolio_return)
d2 <- density(five_CMA$portfolio_return)
plot(range(d1$x, d2$x), range(d1$y, d2$y), type = "n", xlab = "x",
     ylab = "Density")
lines(d1, col = "red")
lines(d2, col = "blue")

plot(five_HML$weights[,2], type = 'l')
