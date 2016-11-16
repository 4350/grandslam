
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
load('data/derived/copula/full_constant_filtered.RData')
load('data/derived/copula/full_dynamic_filtered.RData')

set.seed(1675)
N_RANDOM <- 250000

source('source/copula_simulate.R')

# Simulate stdresid from constant copula ---------------------------------

simulate_stdresid <- function(name, copula) {
  tic(sprintf('stdresid (%s)', name))
  stdresid <- .simulate_stdresid_t(1, copula)
  toc()
  
  colnames(stdresid) <- names(model.GARCH)
  
  filename <- sprintf('data/derived/stdresid/%s.RData', name)
  save(stdresid, file = filename)
}

cl <- makeCluster(spec = detectCores() - 1)
clusterEvalQ(cl, library(devtools))
clusterEvalQ(cl, load_all('australian'))
clusterEvalQ(cl, load_all('wimbledon'))
registerDoParallel(cl)

simulate_stdresid('full_constant_norm', constant_copula_filtered$norm)
simulate_stdresid('full_constant_std', constant_copula_filtered$std)
simulate_stdresid('full_constant_ghst', constant_copula_filtered$ghst)

# Make dynamic copula "constant" to generate unconditional stdresid
dynamic_copula_constant <- lapply(dynamic_copula_filtered, function(filtered) {
  # Make model constant
  filtered$spec@dynamics@alpha <- 0
  filtered$spec@dynamics@beta <- 0
  
  N <- ncol(filtered$spec@dynamics@Omega)
  
  # Null out any dynamic stuff
  filtered$Q <- array(filtered$spec@dynamics@Omega, dim = c(N, N, 1))
  filtered$Correlation <- NULL
  filtered$shocks <- rbind(filtered$shocks[1, ])
  filtered$scores <- NULL
  
  filtered
})

simulate_stdresid('full_dynamic_norm', dynamic_copula_constant$norm)
simulate_stdresid('full_dynamic_std', dynamic_copula_constant$std)
simulate_stdresid('full_dynamic_ghst', dynamic_copula_constant$ghst)
