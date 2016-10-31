#' ARCHIVED BECAUSE OLD COPULA MODELS. 
#' 
#' Simulate 1-step residuals from copula


# Setup ------------------------------------------------------------------

library(ghyp)
library(devtools)
library(tictoc)
load_all('wimbledon')

rm(list = ls())

# First index for which we simulate distribution t + 1
kTStart <- 1
kTEnd <- 2765
kNSim <- 1e4
kModelName <- 'dynamic_ghskt'

output.directory <- file.path('data/derived/stdresid', kModelName)

load(paste0('data/derived/filtered_', kModelName, '.RData'))

dir.create(output.directory, recursive = T, showWarnings = F)

# Test -------------------------------------------------------------------

cluster <- prepare.cluster()

uv <- filtered.series$filtered$copula$uv.dists
garch <- filtered.series$garch

for (t in kTStart:kTEnd) {
  # Shocks for the next period are generated with next period's conditional
  # correlation matrix
  
  tic(paste(kNSim, "simulated for t =", t + 1))
  Correlation_tp1 <- filtered.series$filtered$copula$Correlation[,, t + 1]
  
  mv.dist <- dc.mv.dist(
    filtered.series$copula$params$dist.params,
    Correlation_tp1
  )
  
  stdresid <- shocks2stdresid(ghyp::rghyp(kNSim, mv.dist), uv, garch, cluster)
  write.csv(
    stdresid,
    file.path(
      output.directory,
      paste0(format(kNSim, scientific = F), "_", t + 1, '.csv')
    )
  )
  toc()
}

stopCluster(cluster)
