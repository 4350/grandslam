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
load_all('wimbledon')
rm(list = ls())

MODEL_NAME <- 'dynamic_ghskt'

# Load the distribution and convert to the distributions of simple returns
# Clear the old for memory purposes (the file is approx 1.2 Gb)
load(sprintf('data/derived/distribution_%s.RData', MODEL_NAME))
distribution_simple <- exp(distribution) - 1
rm(distribution)

# Functions --------------------------------------------------------------

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

optimize_cbd <- function(q, distribution_t, eqfun = sum, eqB = 1,
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
                   delta = 1e-6),
    ...
  )
}

# "Unconstrained" Optimization -------------------------------------------

# XXX It may be the case that there is an off-by-one error here; I think
# the first distribution is actually the distribution of t = 2 et cetera.
times <- 1:dim(distribution_simple)[3]
cdb <- rep(NA, length(times))
weights <- matrix(NA, ncol = 6, nrow = length(times))

for (t in times) {
  tic(sprintf('Optimal CDB at t = %d', t))
  op <- optimize_cbd(q = 0.05, distribution_simple[,, t])
  toc()
  
  if (op$convergence > 0) {
    warning(sprintf("No convergence for t = %d", t))
  }
  
  cdb[t] <- -tail(op$values, 1)
  weights[t, ] <- op$pars
}

all_cdb <- cdb
all_weights <- weights

cdb_results <- list(cdb = cdb, weights = weights)
save(cdb_results, file = sprintf('data/derived/cdb_%s_all.RData', MODEL_NAME))

# Only Modern Value -----------------------------------------------------

# Drop classic value from the returns. Note: Replace global distribution
# because it takes so much SPACE!!!
distribution_simple <- distribution_simple[, c(1, 3:6), ]

times <- 1:dim(distribution_simple)[3]
cdb <- rep(NA, length(times))
weights <- matrix(NA, ncol = ncol(distribution_simple), nrow = length(times))

for (t in times) {
  tic(sprintf('Optimal CDB at t = %d', t))
  op <- optimize_cbd(q = 0.05, distribution_simple[,, t])
  toc()
  
  cdb[t] <- -tail(op$values, 1)
  weights[t, ] <- op$pars
}

modern_cdb <- cdb
modern_weights <- weights

cdb_results <- list(cdb = cdb, weights = weights)
save(cdb_results, file = sprintf('data/derived/cdb_%s_modern.RData', MODEL_NAME))

# Only Classic Value ------------------------------------------------------

# Drop RMW and CMA from the distributions. This is a lot faster than having
# constraints in the optimizer
distribution_simple <- distribution_simple[, c(1:4), ]

times <- 1:dim(distribution_simple)[3]
cdb <- rep(NA, length(times))
weights <- matrix(NA, ncol = ncol(distribution_simple), nrow = length(times))

for (t in times) {
  tic(sprintf('Optimal CDB at t = %d', t))
  op <- optimize_cbd(q = 0.05, distribution_simple[,, t])
  toc()
  
  cdb[t] <- -tail(op$values, 1)
  weights[t, ] <- op$pars
}

classic_cdb <- cdb
classic_weights <- weights

cdb_results <- list(cdb = cdb, weights = weights)
save(cdb_results, file = sprintf('data/derived/cdb_%s_classic.RData', MODEL_NAME))


# RMW --------------------------------------------------------------------

distribution_simple <- distribution_simple[, c(1:5), ]

times <- 1:dim(distribution_simple)[3]
cdb <- rep(NA, length(times))
weights <- matrix(NA, ncol = ncol(distribution_simple), nrow = length(times))

for (t in times) {
  tic(sprintf('Optimal CDB at t = %d', t))
  op <- optimize_cbd(q = 0.05, distribution_simple[,, t])
  toc()
  
  cdb[t] <- -tail(op$values, 1)
  weights[t, ] <- op$pars
}

classic_cdb <- cdb
classic_weights <- weights

cdb_results <- list(cdb = cdb, weights = weights)
save(cdb_results, file = sprintf('data/derived/cdb_%s_RMW.RData', MODEL_NAME))

# CMA --------------------------------------------------------------------

distribution_simple <- distribution_simple[, c(1:4, 6), ]

times <- 1:dim(distribution_simple)[3]
cdb <- rep(NA, length(times))
weights <- matrix(NA, ncol = ncol(distribution_simple), nrow = length(times))

for (t in times) {
  tic(sprintf('Optimal CDB at t = %d', t))
  op <- optimize_cbd(q = 0.05, distribution_simple[,, t])
  toc()
  
  cdb[t] <- -tail(op$values, 1)
  weights[t, ] <- op$pars
}

classic_cdb <- cdb
classic_weights <- weights

cdb_results <- list(cdb = cdb, weights = weights)
save(cdb_results, file = sprintf('data/derived/cdb_%s_CMA.RData', MODEL_NAME))


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
do_ew_cdb('dynamic_ghskt', 'classic', c(1:4))
do_ew_cdb('dynamic_ghskt', 'modern', c(1, 3:6))
do_ew_cdb('dynamic_ghskt', 'RMW', c(1:5))
do_ew_cdb('dynamic_ghskt', 'CMA', c(1:4, 6))

# Plotting ---------------------------------------------------------------

rm(list = ls())

MODEL_NAME <- 'dynamic_ghskt'

path_template <- 'data/derived/cdb_%s_%s.RData'
load('data/derived/weekly-full.RData')

load(sprintf(path_template, MODEL_NAME, 'all'))
results_all <- cdb_results

load(sprintf(path_template, MODEL_NAME, 'classic'))
results_classic <- cdb_results

load(sprintf(path_template, MODEL_NAME, 'modern'))
results_modern <- cdb_results

cdb <- data.frame(Week = df$Date[-length(df$Date)],
                  All = results_all$cdb,
                  Classic = results_classic$cdb,
                  Modern = results_modern$cdb)

cdb_long <- gather(cdb, Strategy, CDB, -Week)

g <- ggplot(cdb_long, aes(Week, CDB, colour = Strategy)) +
  geom_line() +
  xlab('')+
  theme_Publication() +
  scale_colour_Publication() +
  coord_cartesian(ylim = c(0.40, 1.00))

OUTPATH <- 'output/CDB/CDB_%s.png'
ggsave(sprintf(OUTPATH, MODEL_NAME), 
       g, device = 'png', width = 14, height = 8, units = 'cm'
)
