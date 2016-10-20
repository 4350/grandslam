# Setup ------------------------------------------------------------------

library(devtools)
library(foreach)
library(doParallel)
load_all('australian')

registerDoParallel(cores = 7)

rm(list = ls())

load('data/derived/weekly-estim.RData')
date <- df.estim$Date
df <- data.frame(df.estim[, -1])
rm(df.estim)

# GARCH Fitting ----------------------------------------------------------

garch_specs <- list(
  Mkt.RF = garch_specgen(0, 0),
  HML = garch_specgen(1, 1),
  SMB = garch_specgen(1, 1),
  Mom = garch_specgen(1, 0),
  RMW = garch_specgen(1, 1),
  CMA = garch_specgen(1, 1)
)

fitted_garch <- garch_fit(garch_specs, df)
stdresid <- foreach(i = seq_along(fitted_garch), .combine = 'cbind') %do% {
  fitted_garch[[i]]@fit$z
}
u <- garch_stdresid2uniform(fitted_garch, stdresid)

# Copula Estimation ------------------------------------------------------

copula <- CopulaSpecification(
  distribution = CopulaDistribution(nu = 8, gamma = rep(0, 6)),
  dynamics = CopulaDynamics(alpha = 0.06, beta = 0.91, phi = 0.10, theta = cbind(rep(0, 6)))
)

# Estimate with time trend
X <- array(
  # Take the time index (seq(nrow(u))), replicate it for each factor and then
  # chop up the matrix into a time-wise array
  matrix(seq(nrow(u)), nrow = ncol(u), ncol = nrow(u), byrow = T),
  c(ncol(u), 1, nrow(u))
)

fitted_copula <- copula_fit(copula, u, X = X,
                            distribution = 't', constant = FALSE,
                            upsilon = TRUE)
