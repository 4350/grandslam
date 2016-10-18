

# Library and setup -------------------------------------------------------
rm(list = ls())
library(tictoc)
library(Rsolnp)
MODEL_NAME = 'dynamic_ghskt'

# Functions ---------------------------------------------------------------

#' Note that this does not get 'standard' tangency, but tangency given that 
#'  weights sum to one and that all weights are positive. 'Standard' tangency
#'  can have extreme weights.

do_optimize_mv <- function(model_name, strategy, selectors) {
  load(sprintf('data/derived/distribution_%s.RData', model_name))
  distribution_simple <- exp(distribution) - 1
  rm(distribution)
  
  # Restrict to the given subset
  colnames(distribution_simple) <- c('Mkt.RF', 'HML', 'SMB', 'Mom', 'RMW', 'CMA')
  distribution_simple <- distribution_simple[, selectors, ]
  
  mv_results <- optimize_mv(distribution_simple)
  rm(distribution_simple)
  
  save(mv_results,
       file = sprintf('data/derived/mv_weights_%s_%s.RData', strategy, model_name))
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


# Get optimal weights for different asset universes -----------------------

do_optimize_mv(MODEL_NAME, 'All', c('Mkt.RF','HML','SMB','Mom','RMW','CMA'))
do_optimize_mv(MODEL_NAME, 'Four+HML', c('Mkt.RF','SMB','Mom','RMW','CMA'))
do_optimize_mv(MODEL_NAME, 'Four+CMA', c('Mkt.RF', 'HML','SMB','Mom','RMW'))

