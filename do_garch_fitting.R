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
#' depending on the model (gjrGARCH or sGARCH) and dist (norm,
#' std, ghst)
#' 
#' @param df.estim 
#' @param model either gjrGARCH or sGARCH
#' @param dist either norm, std or ghst
#' 
#' @return saves the garch_fits objects in derived data
do_garch_fit <- function(df.estim, model, dist, nOOS = 0) {
  
  # Estimate GARCH specifications ----
  speclist <- list(
    ARMA00 = garch.specgen(0,0, model = model, dist = dist),
    ARMA10 = garch.specgen(1,0, model = model, dist = dist),
    ARMA01 = garch.specgen(0,1, model = model, dist = dist),
    ARMA11 = garch.specgen(1,1, model = model, dist = dist),
    ARMA20 = garch.specgen(2,0, model = model, dist = dist),
    ARMA21 = garch.specgen(2,1, model = model, dist = dist),
    ARMA22 = garch.specgen(2,2, model = model, dist = dist),
    ARMA02 = garch.specgen(0,2, model = model, dist = dist),
    ARMA12 = garch.specgen(1,2, model = model, dist = dist),
    ARMA30 = garch.specgen(3,0, model = model, dist = dist),
    ARMA31 = garch.specgen(3,1, model = model, dist = dist),
    ARMA32 = garch.specgen(3,2, model = model, dist = dist),
    ARMA33 = garch.specgen(3,3, model = model, dist = dist),
    ARMA03 = garch.specgen(0,3, model = model, dist = dist),
    ARMA13 = garch.specgen(1,3, model = model, dist = dist),
    ARMA23 = garch.specgen(2,3, model = model, dist = dist)
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
  if(nOOS != 0) {
    save(fits, file = sprintf('data/derived/garch/oos_garch_fits_%s_%s.RData', model, dist))
  } else {
    save(fits, file = sprintf('data/derived/garch/garch_fits_%s_%s.RData', model, dist))
  }
}

do_garch_fit(df.estim, 'gjrGARCH', 'ghst')
do_garch_fit(df.estim, 'gjrGARCH', 'std')
do_garch_fit(df.estim, 'gjrGARCH', 'norm')

do_garch_fit(df.estim, 'sGARCH', 'ghst')
do_garch_fit(df.estim, 'sGARCH', 'std')
do_garch_fit(df.estim, 'sGARCH', 'norm')


# Garch OOS ---------------------------------------------------------------

do_garch_fit(df.estim, 'gjrGARCH', 'ghst', nOOS = 914) #first date is 1999-01-01, approx 1/3 of full sample
do_garch_fit(df.estim, 'sGARCH', 'ghst', nOOS = 914)