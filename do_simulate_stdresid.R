
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
load('data/derived/copula/legacy/oos_constant_filtered.RData')

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

registerDoParallel(cores = 7)
simulate_stdresid('full_constant_norm', constant_copula_filtered$norm)
simulate_stdresid('full_constant_std', constant_copula_filtered$std)
simulate_stdresid('full_constant_ghst', constant_copula_filtered$ghst)
