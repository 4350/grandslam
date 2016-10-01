#' Get the sigma matrix over time for tail dependency use
#' takes hard coded parameters
#' 

# Reset ----

rm(list = ls())
load('data/derived/garch_unires_model.RData')

# Split u into modern and classic
u <- df.u[, -1]
rm(df.u)

u.modern <- u[,c('Mkt.RF','SMB','Mom','RMW','CMA')]
u.classic <- u[,c('Mkt.RF','HML','SMB','Mom')]

# Libraries ----
library(ghyp)
library(parallel)
library(devtools)
load_all('wimbledon')

# Parameter load
# Modern
load('data/derived/dcopula_param_ghskt_modern.RData')
dist.params.modern <- list(
  df = param.ghskt$par[1],
  skew = param.ghskt$par[2:(length(param.ghskt$par)-2)]
)

alpha.modern = tail(param.ghskt$par,2)[1]
beta.modern = tail(param.ghskt$par,2)[2]

# Classic
load('data/derived/dcopula_param_ghskt_classic.RData')
dist.params.classic <- list(
  df = param.ghskt$par[1],
  skew = param.ghskt$par[2:(length(param.ghskt$par)-2)]
)

alpha.classic = tail(param.ghskt$par,2)[1]
beta.classic = tail(param.ghskt$par,2)[2]

rm(param.ghskt)
rm(u)

# Initiate clusters
prepare.cluster <- function() {
  cluster <- makeCluster(detectCores() - 1)
  clusterEvalQ(cluster, library(ghyp))
  clusterEvalQ(cluster, library(devtools))
  clusterEvalQ(cluster, load_all('wimbledon'))
  cluster
}

cluster <- prepare.cluster()

get.sigma <- function(u, dist.params, alpha, beta, cluster) {
  # Build univariate distributions as they're used for the construction
  # of our shocks; the MV distribution is built for each t based on the
  # Correlation matrix generated
  uv.dists <- lapply(seq(ncol(u)), function(i) {
    # Gaussian Case
    if (is.null(dist.params$df)) {
      dist <- ghyp::gauss()
    }
    # Symmetric T case
    else if (is.null(dist.params$skew)) {
      dist <- student.t(nu = dist.params$df)
    }
    # Skewed T case
    else {
      dist <- student.t(nu = dist.params$df, gamma = dist.params$skew[i])
    }
  })
  shocks <- dc.shocks(u, uv.dists, cluster)
  
  # Compute the DCC model parameters
  shocks.std <- dc.shocks.std(shocks, uv.dists, alpha, beta)
  Omega <- dc.Omega(shocks.std)
  Q <- dc.Q(shocks.std, Omega, alpha, beta)
  Correlation <- dc.Correlation(Q, cluster = cluster)
  
}

Correlations.modern <- get.sigma(u.modern, dist.params.modern, alpha.modern, beta.modern, cluster = cluster)
Correlations.classic <- get.sigma(u.classic, dist.params.classic, alpha.classic, beta.classic, cluster = cluster)
  
# Use correlation matrices and shocks of arbitrary small u to find tail dependence over time



pghyp(q = c(0.01))