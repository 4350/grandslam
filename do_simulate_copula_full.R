# Setup ------------------------------------------------------------------

rm(list = ls())

library(foreach)
library(doParallel)
library(tictoc)
library(devtools)
library(rugarch)
library(memoise)
load_all('australian')
load_all('wimbledon')

load('data/derived/garch/model_GARCH_chosen.RData')
load('data/derived/weekly-estim.RData')
load('data/derived/garch/model_GARCH_chosen_filtered.RData')

source('source/copula_simulate.R')
set.seed(1675)

# Needed for GARCH path
specs <- garch.fit2spec(model.GARCH)

# Filtered data to rebuild GARCH paths
filtered_series <- list(
  series = data.frame(df.estim[, -1]),
  resid = sapply(filtered, rugarch::residuals),
  sigma = sapply(filtered, rugarch::sigma)
)
rm(filtered)

# Configuration
N_RANDOM = 1e4
T_SEQ <- seq(1, 2765) # Times WHEN distributions are simulated FOR t + 1

# Do Simulation ----------------------------------------------------------
cl <- makeCluster(spec = detectCores() - 1)
clusterEvalQ(cl, library(devtools))
clusterEvalQ(cl, load_all('wimbledon'))
clusterEvalQ(cl, load_all('australian'))
registerDoParallel(cl)
# registerDoParallel(cores = 7)

load('data/derived/copula/full_constant_filtered.RData')
do_simulate('full_constant_norm', constant_copula_filtered$norm)
do_simulate('full_constant_std',  constant_copula_filtered$std)
do_simulate('full_constant_ghst', constant_copula_filtered$ghst)

load('data/derived/copula/full_dynamic_filtered.RData')
do_simulate('full_dynamic_norm', dynamic_copula_filtered$norm)
do_simulate('full_dynamic_std',  dynamic_copula_filtered$std)
do_simulate('full_dynamic_ghst', dynamic_copula_filtered$ghst)

stopCluster()