# Fit a battery of GARCH models to Fama-French factors
#
# Prerequisites
#
# * Need to run `do_load.R` first to properly prepare sample

# Libraries ----
library(rugarch)
library(parallel)
library(devtools)
load_all('wimbledon')

# Reset Workspace ----
rm(list = ls())
load('data/derived/weekly-estim.RData')

# Estimate GARCH specifications ----
speclist <- list(
  ARMA00 = garch.specgen(0,0),
  ARMA10 = garch.specgen(1,0),
  ARMA01 = garch.specgen(0,1),
  ARMA11 = garch.specgen(1,1),
  ARMA20 = garch.specgen(2,0),
  ARMA21 = garch.specgen(2,1),
  ARMA22 = garch.specgen(2,2),
  ARMA02 = garch.specgen(0,2),
  ARMA12 = garch.specgen(1,2),
  ARMA30 = garch.specgen(3,0),
  ARMA31 = garch.specgen(3,1),
  ARMA32 = garch.specgen(3,2),
  ARMA33 = garch.specgen(3,3),
  ARMA03 = garch.specgen(0,3),
  ARMA13 = garch.specgen(1,3),
  ARMA23 = garch.specgen(2,3)
)

# Parallelize the estimation of GARCH, by splitting up the estimation of each
# GARCH in parallel steps.
cluster <- makeCluster(detectCores() - 1)
clusterEvalQ(cluster, library(rugarch))

fits <- lapply(
  df.estim[, -1],
  function(factor) {
    parLapply(
      cluster,
      speclist,
      function(spec, factor) {
        ugarchfit(spec, factor, model = "hybrid")
      },
      factor = factor
    )
  }
)
stopCluster(cluster)
rm(cluster)

save(fits, file = 'data/derived/garch_fits.RData')
  