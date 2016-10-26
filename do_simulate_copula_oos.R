#' Simulate from copula models OOS
#' 
#' These are the distributions of returns given models estimated on a
#' restricted sample.
#' 
#' See `do_simulate_independence.R` for more straightforward code that ignores
#' copula. Before running this code, I need fully filtered data for each
#' copula (which comes after the copula estimation in `do_fit_copula2.R`)


# Setup ------------------------------------------------------------------

library(foreach)
library(doParallel)
library(tictoc)
library(devtools)
library(rugarch)
library(memoise)
load_all('australian')
load_all('wimbledon')

registerDoParallel(cores = 7)

rm(list = ls())

N_RANDOM <- 1e4          # Number of simulation paths per t
T_SEQ <- seq(1852, 2765) # Times WHEN distributions are simulated FOR t + 1

set.seed(1675)

# The filtered data needs to be in a simpler format, we basically expect
# it to be a list of matrices series, sigma and resid
load('data/derived/garch/oos_model_GARCH_chosen_filtered.RData')
load('data/derived/weekly-estim.RData')
filtered_series <- list(
  series = data.frame(df.estim[, -1]),
  resid = sapply(filtered, rugarch::residuals),
  sigma = sapply(filtered, rugarch::sigma)
)
rm(filtered)

# Needed for GARCH path
load('data/derived/garch/oos_model_GARCH_chosen.RData')
specs <- garch.fit2spec(model.GARCH)

# FUNCTIONS

garch_path_t <- function(t, filtered_copula) {
  # Get simulated stdresid for this ("next") period
  stdresid_t <- .simulate_stdresid_t(t, filtered_copula)
  
  foreach (i = seq_along(model.GARCH), .combine = 'cbind') %do% {
    path <- garch.path.it(i, t, specs, filtered_series, stdresid_t)
    t(path@path$seriesSim)
  }
}

#' Simulate standardized residuals from copula
#' 
#' Uses a memoised version of the underlying simulator if using a constant
#' copula.
.simulate_stdresid_t <- function(t, filtered_copula) {
  if (copula_is_constant(filtered_copula$spec)) {
    return(.do_simulate_stdresid_t_constant(filtered_copula))
  }
  
  .do_simulate_stdresid_t(t, filtered_copula)
}

#' Simulate t+1 uniforms from copula
.copula_simulate_t <- function(t, filtered_copula) {
  Q_t <- filtered_copula$Q[,, t]
  shocks_t <- filtered_copula$shocks[t, ]
  
  t(simplify2array(
    copula_simulate(filtered_copula$spec, 1, N_RANDOM, Q_t, shocks_t)
  ))
}

# Do the actual simulation
.do_simulate_stdresid_t <- function(t, filtered_copula) {
  tic('simulate_t')
  u <- .copula_simulate_t(t, filtered_copula)
  toc()
  tic('uniform2stdresid')
  x <- garch_uniform2stdresid(model.GARCH, u)
  toc()
  x
}

# Memoised version for constant 
.do_simulate_stdresid_t_constant <- memoise(
  function(c) .do_simulate_stdresid_t(1, c)
)

do_simulate <- function(name, copula) {
  distribution <- array(NA, dim = c(N_RANDOM, 6, length(T_SEQ)))
  colnames(distribution) <- names(model.GARCH)
  
  for (t_i in seq_along(T_SEQ)) {
    # Select INDEX because we don't want to fill up with all the unsimulated
    # periods
    t <- T_SEQ[t_i]
    
    tic(sprintf('%s: t = %d', name, t))
    distribution[,, t_i] <- garch_path_t(t, copula)
    toc()
  }
  
  file <- file.path(
    'data/derived/distributions',
    sprintf('oos_%s_%d.RData', name, N_RANDOM)
  )
  
  save(distribution, file = file)
  
  # avoid returning the results of save in case to allow for garbage collection
  return('OK')
}

# Run Simulations --------------------------------------------------------

load('data/derived/copula/oos_constant_filtered.RData')
do_simulate('constant_norm', constant_copula_filtered$norm)
do_simulate('constant_std',  constant_copula_filtered$std)
do_simulate('constant_ghst', constant_copula_filtered$ghst)

load('data/derived/copula/oos_dynamic_filtered.RData')
do_simulate('dynamic_norm', dynamic_copula_filtered$norm)
do_simulate('dynamic_std',  dynamic_copula_filtered$std)
do_simulate('dynamic_ghst', dynamic_copula_filtered$ghst)
