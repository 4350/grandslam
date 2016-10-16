#' We compute the CDB measure of Christoffersen (2012) using distributions
#' of asset returns. Note that the distributions themselves are for log-returns
#' but we compute ES/CDB based on simple returns.
#' 
#' This is a minor thing in the grander scheme of things -- but it does ensure
#' the weighted average property holds.

# Setup ------------------------------------------------------------------

library(tictoc)
library(Rsolnp)
library(ggplot2)
library(tidyr)
library(devtools)
library(foreach)
load_all('wimbledon')
rm(list = ls())

MODEL_NAME <- 'dynamic_ghskt'

# Functions 

#' Compute VaR, ES and CDB for a portfolio
#'
#' @param q quantile
#' @param weights N weights in assets
#' @param returns MxN distribution of returns matrix
#'
#' @return List with portfolio VaR (var), ES (es) and CDB (cdb)
#' @export
risk_measures <- function(weights, q, returns) {
  # Note: This does not depend on weights and could therefore be computed once
  # so that we're only optimizing for the portfolio. It doesn't take long,
  # however, only about 7 milliseconds with 6 * 1e5 returns
  vars <- apply(returns, 2, function(r) -quantile(r, q))
  es <- lapply(seq_along(vars), function(i) {
    r <- returns[, i]
    -mean(r[r <= -vars[i]])
  })
  es <- unlist(es)
  
  portfolio_returns <- returns %*% weights
  portfolio_var <- -quantile(portfolio_returns, q, names = FALSE)
  portfolio_es <- -mean(portfolio_returns[portfolio_returns <= -portfolio_var])
  
  es_ub <- sum(weights * es)
  es_lb <- portfolio_var
  
  cdb <- (es_ub - portfolio_es) / (es_ub - es_lb)
  
  list(
    var = portfolio_var,
    es = portfolio_es,
    cdb = cdb
  )
}

optimize_cdb <- function(q, distribution_t, eqfun = sum, eqB = 1,
                         LB = NULL, x0 = NULL, ...) {
  N <- ncol(distribution_t)
  
  if (is.null(LB)) {
    LB <- rep(0, N)
  }
  
  # Start with EW portfolio
  if (is.null(x0)) {
    x0 <- rep(1 / N, N)
  }
  
  fn <- function(weights) -risk_measures(weights, q, distribution_t)$cdb
  
  solnp(
    x0,
    fn,
    eqfun = eqfun,
    eqB = eqB,
    LB = rep(0, N),
    
    # Optimizer keeps climbing the multidimensional mountains with tiny
    # tiny steps. Take some strides (default delta = 1e-7)!!!
    control = list(trace = 0,
                   delta = 1e-5),
    ...
  )
}

best_cdb <- function(distribution) {
  times <- 1:dim(distribution)[3]
  cdb <- rep(NA, length(times))
  weights <- matrix(NA, ncol = ncol(distribution), nrow = length(times))
  
  for (t in times) {
    tic(sprintf('Optimal CDB at t = %d', t))
    op <- optimize_cdb(q = 0.05, distribution[,, t])
    toc()
    
    if (op$convergence > 0) {
      warning(sprintf('No convergence for t = %d', t))
    }
    
    cdb[t] <- -tail(op$values, 1)
    weights[t, ] <- op$pars
  }
  
  list(
    cdb = cdb,
    weights = weights
  )
}

do_best_cdb <- function(model_name, strategy, selectors) {
  load(sprintf('data/derived/distribution_%s.RData', model_name))
  distribution_simple <- exp(distribution) - 1
  rm(distribution)
  
  # Restrict to the given subset
  colnames(distribution_simple) <- c('Mkt.RF', 'HML', 'SMB', 'Mom', 'RMW', 'CMA')
  distribution_simple <- distribution_simple[, selectors, ]
  
  cdb_results <- best_cdb(distribution_simple)
  rm(distribution_simple)
  
  save(cdb_results,
       file = sprintf('data/derived/cdb_%s_%s.RData', strategy, model_name))
}

# CDB Optimization -------------------------------------------------------

do_best_cdb(MODEL_NAME, 'all',
            c('Mkt.RF', 'HML', 'SMB', 'Mom', 'RMW', 'CMA'))

do_best_cdb(MODEL_NAME, 'modern',
            c('Mkt.RF', 'SMB', 'Mom', 'RMW', 'CMA'))

do_best_cdb(MODEL_NAME, 'HML', c('Mkt.RF', 'HML', 'SMB', 'Mom'))
do_best_cdb(MODEL_NAME, 'RMW', c('Mkt.RF', 'SMB', 'Mom', 'RMW'))
do_best_cdb(MODEL_NAME, 'CMA', c('Mkt.RF', 'SMB', 'Mom', 'CMA'))

do_best_cdb(MODEL_NAME, 'RMW+HML', c('Mkt.RF', 'SMB', 'Mom', 'RMW', 'HML'))
do_best_cdb(MODEL_NAME, 'RMW+CMA', c('Mkt.RF', 'SMB', 'Mom', 'RMW', 'CMA'))

# Equal Weights ----------------------------------------------------------

ew_cdb <- function(distribution) {
  times <- 1:dim(distribution)[3]
  weights <- rep(1 / ncol(distribution), ncol(distribution))
  
  cdb <- rep(NA, length(times))
  
  foreach(t = times, .combine = 'c') %do% {
    risk_measures(weights = weights,
                  q = 0.05,
                  returns = distribution[,, t])$cdb
  }
}

do_ew_cdb <- function(model_name, name, factors) {
  load(sprintf('data/derived/distribution_%s.RData', model_name))
  distribution <- distribution[, factors, ]
  
  distribution_simple <- exp(distribution) - 1
  rm(distribution)
  
  cdb <- ew_cdb(distribution_simple)
  rm(distribution_simple)
  
  filename <- sprintf('data/derived/cdb_%s_%s_ew.RData', model_name, name)
  save(cdb, file = filename)
}

do_ew_cdb('dynamic_ghskt', 'all', c(1:6))
do_ew_cdb('dynamic_ghskt', 'HML', c(1:4))
do_ew_cdb('dynamic_ghskt', 'modern', c(1, 3, 4, 5, 6))
do_ew_cdb('dynamic_ghskt', 'RMW', c(1, 3, 4, 5))
do_ew_cdb('dynamic_ghskt', 'CMA', c(1, 3, 4, 6))

do_ew_cdb('dynamic_ghskt', 'RMW+HML', c(1, 2, 3, 4, 5   ))
do_ew_cdb('dynamic_ghskt', 'RMW+CMA', c(1,    3, 4, 5, 6))

# Plotting ---------------------------------------------------------------

rm(list = ls())

MODEL_NAME <- 'dynamic_ghskt'

load('data/derived/weekly-full.RData')

load_cdb_optim <- function(name) {
  load(sprintf('data/derived/cdb_%s_%s.RData', MODEL_NAME, name))
  cdb_results$cdb
}

load_cdb_ew <- function(name) {
  load(sprintf('data/derived/cdb_%s_%s_ew.RData', MODEL_NAME, name))
  cdb
}

load_cdb <- function(name) {
  optim <- load_cdb_optim(name)
  ew <- load_cdb_ew(name)
  
  data.frame(
    Week = df$Date[-length(df$Date)],
    Optimized = load_cdb_optim(name),
    EW = load_cdb_ew(name)
  )
}

plot_cdb <- function(cdb, name, width, height) {
  cdb_long <- gather(cdb, Routine, CDB, Optimized, factor_key = TRUE)
  
  g <- ggplot(cdb_long, aes(Week, CDB, colour = Strategy)) +
    geom_line() +
    xlab('')+
    theme_Publication() +
    scale_colour_Publication() +
    coord_cartesian(ylim = c(0.40, 1.00))
  
  
  ggsave(filename = sprintf('output/CDB/cdb--%s.png', name),
         plot = g, units = c('cm'), width = width, height = height)
  g
}

# All/Modern/Classic -----------------------------------------------------

cdb <- bind_rows(
  All = load_cdb('all'),
  Modern = load_cdb('modern'),
  Classic = load_cdb('HML'),
  .id = 'Strategy'
)

plot_cdb(cdb, 'modern-classic', width = 14, height = 7)

# RMW + CMA or HML -------------------------------------------------------

cdb <- bind_rows(
  `RMW + CMA` = load_cdb('RMW+CMA'),
  `RMW + HML` = load_cdb('RMW+HML'),
  .id = 'Strategy'
)

plot_cdb(cdb, 'rmw_cma-rmw_hml', width = 14, height = 7)


# 3-factors --------------------------------------------------------------

cdb <- bind_rows(
  RMW = load_cdb('RMW'),
  CMA = load_cdb('CMA'),
  HML = load_cdb('HML'),
  .id = 'Strategy'
)

plot_cdb(cdb, 'rmw-cma-hml', width = 14, height = 7)
