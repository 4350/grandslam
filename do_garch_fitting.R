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

#' Function to do GARCH fit and save to different objects
#' depending on the submodel (if GJR) and dist (norm,
#' std, ghst)
#' 
#' @param df.estim 
#' @param submodel either GJRGARCH or GARCH
#' @param dist either norm, std or ghst
#' 
#' @return saves the garch_fits objects in derived data
do_garch_fit <- function(df.estim, submodel, dist){
  # Set OOS to zero
  nOOS <- 0
  
  # Estimate GARCH specifications ----
  speclist <- list(
    ARMA00 = garch.specgen(0,0, submodel = submodel, dist = dist),
    ARMA10 = garch.specgen(1,0, submodel = submodel, dist = dist),
    ARMA01 = garch.specgen(0,1, submodel = submodel, dist = dist),
    ARMA11 = garch.specgen(1,1, submodel = submodel, dist = dist),
    ARMA20 = garch.specgen(2,0, submodel = submodel, dist = dist),
    ARMA21 = garch.specgen(2,1, submodel = submodel, dist = dist),
    ARMA22 = garch.specgen(2,2, submodel = submodel, dist = dist),
    ARMA02 = garch.specgen(0,2, submodel = submodel, dist = dist),
    ARMA12 = garch.specgen(1,2, submodel = submodel, dist = dist),
    ARMA30 = garch.specgen(3,0, submodel = submodel, dist = dist),
    ARMA31 = garch.specgen(3,1, submodel = submodel, dist = dist),
    ARMA32 = garch.specgen(3,2, submodel = submodel, dist = dist),
    ARMA33 = garch.specgen(3,3, submodel = submodel, dist = dist),
    ARMA03 = garch.specgen(0,3, submodel = submodel, dist = dist),
    ARMA13 = garch.specgen(1,3, submodel = submodel, dist = dist),
    ARMA23 = garch.specgen(2,3, submodel = submodel, dist = dist)
  )
  
  # Parallelize the estimation of GARCH, by splitting up the estimation of each
  # GARCH in parallel steps.
  cluster <- makeCluster(detectCores() - 1)
  clusterEvalQ(cluster, library(rugarch))
  
  fits <- lapply(
    df.estim[, -1],
    function(factor, nOOS) {
      parLapply(
        cluster,
        speclist,
        function(spec, factor, nOOS) {
          ugarchfit(spec, factor, out.sample = nOOS, model = "hybrid")
        },
        factor = factor,
        nOOS = nOOS
      )
    },
    nOOS = nOOS
  )
  stopCluster(cluster)
  rm(cluster)  
  save(fits, file = sprintf('data/derived/garch_fits_%s_%s.RData', submodel, dist))
}

do_garch_fit(df.estim, 'GJRGARCH', 'ghst')
do_garch_fit(df.estim, 'GJRGARCH', 'std')
do_garch_fit(df.estim, 'GJRGARCH', 'norm')

do_garch_fit(df.estim, 'GARCH', 'ghst')
do_garch_fit(df.estim, 'GARCH', 'std')
do_garch_fit(df.estim, 'GARCH', 'norm')
