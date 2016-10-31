#' ARCHIVED BECAUSE OLD COPULA MODELS

#' One-step ahead distributions
#' 
#' Given information available at time t, we can simulate the distribution of
#' returns for t+1.


# Setup ------------------------------------------------------------------

library(devtools)
library(rugarch)
library(parallel)
library(tictoc)
load_all('wimbledon')

rm(list = ls())

kModelName <- 'dynamic_ghskt'
kNSim <- 1e4

# Load filtered dataset (and models)
load(paste0('data/derived/filtered_', kModelName, '.RData'))

# Simulate 1-step distributions ------------------------------------------

garch.path.t <- function(t) {
  # Load standardized residuals for this model/time from file
  csv.file <- sprintf("%d_%d.csv", kNSim, t)
  stdresid <- read.csv(
    file.path('data/derived/stdresid', kModelName, csv.file),
    row.names = 1
  )
  
  # Use standardized residuals to generate 1-step distribution of series
  sapply(seq(ncol(stdresid)), function(i) {
    path <- garch.path.it(
      i,
      t,
      filtered.series$garch,
      filtered.series$filtered,
      stdresid
    )
    
    # We only care about the series path; the sigma/residuals are irrelevant
    # and the rseed and model to
    t(path@path$seriesSim)
  })
}


# Run for each t ---------------------------------------------------------

tt <- 2:2766
distribution <- array(
  NA,
  dim = c(
    kNSim,
    length(filtered.series$garch),
    length(tt)
  )
)

for (ti in seq(tt)) {
  t <- tt[ti]
  tic(paste("Distribution for t =", t))
  distribution[,, ti] <- garch.path.t(t)
  toc()
}

save(
  distribution,
  file = file.path(
    'data/derived',
    sprintf('distribution_%s.RData', kModelName)
  )
)