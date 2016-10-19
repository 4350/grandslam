

# Library and setup -------------------------------------------------------
rm(list = ls())
library(tictoc)
library(Rsolnp)
MODEL_NAME = 'dynamic_ghskt'

# Functions ---------------------------------------------------------------

#' Note that this does not get 'standard' tangency, but tangency given that 
#'  weights sum to one and that all weights are positive. 'Standard' tangency
#'  can have extreme weights.

do_optimize_mv <- function(model_name, selectors) {
  # Load distribution simulated data
  load(sprintf('data/derived/distribution_%s.RData', model_name))
  distribution_simple <- exp(distribution) - 1
  rm(distribution)
  
  # Restrict to the given subset
  colnames(distribution_simple) <- c('Mkt.RF', 'HML', 'SMB', 'Mom', 'RMW', 'CMA')
  distribution_simple <- distribution_simple[, selectors, ]
  
  mv_results <- optimize_mv(distribution_simple)
  rm(distribution_simple)
  
  return(mv_results)
}

# Give me distribution 3 dimensional array and i give you the mv weights
optimize_mv <- function(distribution) {
  times <- 1:dim(distribution)[3]
  N <- ncol(distribution)
  T <- length(times)
  weights <- matrix(NA, ncol = N, nrow = T)
  sr <- rep(NA, T)
  
  for (t in times) {
    tic(sprintf('Optimal weights at t = %d', t))
    op <- optimize_mv_1_period(distribution[,,t])
    toc()
    
    if (op$convergence > 0) {
      warning(sprintf('No convergence for t = %d', t))
    }
    
    weights[t, ] <- op$pars
    sr[t] <- -tail(op$values,1)
  }
  
  # Calculate portfolio returns
  list(
    weights = weights, sr = sr
  )
}

# Returns optimized objet from solnp
optimize_mv_1_period <- function(distribution_t, eqfun = sum, eqB = 1,
                                 LB = NULL, x0 = NULL, ...) {
  
  N <- ncol(distribution_t)
  # No negative weights
  if (is.null(LB)) {
    LB <- rep(0, N)
  }
  
  # Start with EW portfolio
  if (is.null(x0)) {
    x0 <- rep(1 / N, N)
  }
  
  # Inputs needed from simulated 1-step-ahead distribution
  mu <- colMeans(distribution_t)
  sigma <- cov(distribution_t)
  
  fn <- function(weights, mu, sigma) -sharpe_ratio(weights, mu, sigma)[1]
  
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
    mu = mu,
    sigma = sigma
  )
}

sharpe_ratio <- function(weights, mu, sigma) {
  # 
  sr <- (weights %*% mu) / sqrt(weights %*% sigma %*% weights)
}


get_mv_portfolio_results <- function(model_name, strategy, selectors, realized) {
  # Load MV results and SR
  mv_results <- do_optimize_mv(model_name, selectors)
  colnames(mv_results$weights) <- selectors
  
  T = nrow(mv_results$weights)
  N = ncol(mv_results$weights)
  
  # Get realized returns and dates
  load('data/derived/weekly-estim.RData')
  dates <- tail(df.estim[,'Date'], T)
  realized <- tail(df.estim[, selectors], T)
  rm(df.estim)
  
  # Calculate portfolio return
  portfolio_return <- rowSums(mv_results$weights * realized)
  
  out.list <- list(
    Date = dates,
    sr = mv_results$sr,
    weights = mv_results$weights,
    portfolio_return = portfolio_return
  )
  
  save(out.list, file = sprintf('data/derived/mv_results_%s_%s.Rdata', strategy, model_name))
}




# Get results for different portfolios ------------------------------------

get_mv_portfolio_results(MODEL_NAME, 'All', c('Mkt.RF','HML','SMB','Mom','RMW','CMA'))
get_mv_portfolio_results(MODEL_NAME, 'Four+HML', c('Mkt.RF', 'HML', 'SMB', 'Mom', 'RMW'))
get_mv_portfolio_results(MODEL_NAME, 'Four+CMA', c('Mkt.RF', 'SMB', 'Mom', 'RMW', 'CMA'))


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

m(five_CMA$weights[,5], type = 'l')
plot(five_HML$weights[,2], type = 'l')
