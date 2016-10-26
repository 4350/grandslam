#' Independence Copula Simulation
#' 
#' If shocks are completely unrelated, then the source of randomness is
#' really just a cloud of random uniforms, which get upjected into GARCH
#' and run through the model.
#' 

# Setup ------------------------------------------------------------------

library(foreach)
library(tictoc)
library(devtools)
library(rugarch)
load_all('wimbledon')

rm(list = ls())

# Number of simulation paths to do
N_RANDOM <- 1e5

# First time for which we simulate t + 1
T_SEQ <- seq(1852, 2765)

# Use the same set of runif each period
set.seed(403)
u <- foreach(i = seq(6), .combine = 'cbind') %do% runif(N_RANDOM)

# Find the shocks (same distribution every period)
load('data/derived/garch/oos_model_GARCH_chosen.RData')

stdresid <- foreach(i = seq(6), .combine = 'cbind') %do% {
  garch.qghyp.rugarch(
    u[, i],
    shape = model.GARCH[[i]]@fit$coef['shape'],
    skew = model.GARCH[[i]]@fit$coef['skew']
  )
}

rm(u, i)

# The filtered data needs to be in a simpler format, we basically expect
# it to be a list of matrices series, sigma and resid
load('data/derived/garch/oos_model_GARCH_chosen_filtered.RData')
load('data/derived/weekly-estim.RData')
filtered <- list(
  series = data.frame(df.estim[, -1]),
  resid = sapply(filtered, rugarch::residuals),
  sigma = sapply(filtered, rugarch::sigma)
)

# Needed for GARCH path
specs <- garch.fit2spec(model.GARCH)

# Simulation -------------------------------------------------------------

#' Run 1-step GARCH path for each GARCH for each t
#' 
#' Note: Uses the *same* stdresid for each run!!!
garch.path.t <- function(t) {
  foreach(i = seq_along(model.GARCH), .combine = 'cbind') %do% {
    path <- garch.path.it(i, t, specs, filtered, stdresid)
    t(path@path$seriesSim)
  }
}

# Do for each t
distribution <- array(NA, dim = c(N_RANDOM, 6, length(T_SEQ)))
colnames(distribution) <- names(model.GARCH)

for (t_i in seq_along(T_SEQ)) {
  # Select INDEX because we don't want to fill up with all the unsimulated
  # periods
  t <- T_SEQ[t_i]
  
  tic(sprintf('Distribution t = %d', t))
  distribution[,, t_i] <- garch.path.t(t)
  toc()
}

save(
  distribution,
  file = file.path(
    'data/derived/distributions',
    sprintf('oos_indep_%d.RData', N_RANDOM)
  )
)
