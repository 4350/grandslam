# Perform grid search optimization of dynamic copula model

# Load Libraries ----
library(parallel)
library(compiler)
library(ghyp)
loadcmp('func/dynamicCopula.Rc')

# Reset Environment to GARCH uniform residuals ----
rm(list = ls())
load('data/derived/garch5f-u.RData')

# Log Likelihood Functions ----

grid.optimJoint <- function(u, df, skew, cluster) {
  # Function to be optimized
  objectiveFn <- function(params, shocks, uv.dists, mv.df, mv.skew, marginalLL, cluster = cluster) {
    alpha = params[1]
    beta = params[2]
    
    jointLL <- jointLogLikelihood(
      shocks,
      uv.dists,
      mv.df = mv.df,
      mv.skew = mv.skew,
      alpha = alpha,
      beta = beta,
      cluster = cluster
    )
    
    -(sum(jointLL) - marginalLL)
  }
  
  uv.dists <- lapply(skew, function(sk) student.t(nu = df, gamma = sk))
  shocks <- dc.shocks(u, uv.dists, cluster)
  marginalLL <- sum(rowSums(marginalLogLikelihood(shocks, uv.dists, cluster)))
  
  # Only optimize for params controlling the MV distribution
  params <- c(0.03, 0.96)
  
  # objectiveFn(params, shocks, uv.dists, df, skew, marginalLL, cluster)

  constrOptim(
    params,
    objectiveFn,
    grad = NULL,
    ui = rbind(
      c( 1,  0),
      c( 0,  1),
      c(-1, -1)
    ),
    ci = rbind(
      0,
      0,
      -0.9999
    ),
    control = list(
      trace = 6
    ),

    shocks = shocks,
    uv.dists = uv.dists,
    mv.df = df,
    mv.skew = skew,
    marginalLL = marginalLL,
    cluster = cluster
  )
}

# Test ----
kNumCores <- detectCores() - 1
cluster <- makeCluster(kNumCores)
clusterEvalQ(cluster, library(ghyp))
clusterExport(cluster, "marginalLogLikelihood")
clusterExport(cluster, "jointLogLikelihood")
clusterExport(cluster, "dc.shocks")
clusterExport(cluster, "dc.shocks.std")
clusterExport(cluster, "dc.Q")
clusterExport(cluster, "dc.Omega")
clusterExport(cluster, "dc.Correlation")

ptm <- proc.time()
grid.optimJoint(u, df = 6, skew = rep(0, 6), cluster = cluster)
proc.time() - ptm

stopCluster(cluster)

# Test of total ll ----
kNumCores <- detectCores() - 1
cluster <- makeCluster(kNumCores)
clusterEvalQ(cluster, library(ghyp))
clusterExport(cluster, "marginalLogLikelihood")
clusterExport(cluster, "jointLogLikelihood")
clusterExport(cluster, "dc.shocks")
clusterExport(cluster, "dc.shocks.std")
clusterExport(cluster, "dc.Q")
clusterExport(cluster, "dc.Omega")
clusterExport(cluster, "dc.Correlation")

ptm <- proc.time()
totalLogLikelihood(u, 6, skew = rep(0, 6), 0.03, 0.96, cluster)
proc.time() - ptm

stopCluster(cluster)

# TIMING TEST ----
kNumCores <- detectCores() - 1
cluster <- makeCluster(kNumCores)
clusterEvalQ(cluster, library(ghyp))
clusterExport(cluster, "marginalLogLikelihood")
clusterExport(cluster, "jointLogLikelihood")
clusterExport(cluster, "dc.shocks")
clusterExport(cluster, "dc.shocks.std")
clusterExport(cluster, "dc.Q")
clusterExport(cluster, "dc.Omega")
clusterExport(cluster, "dc.Correlation")


df = 6
skew = rep(0, 6)
alpha = 0.03
beta = 0.96

# GETTING SHOCKS TAKES HOW LONG?
"SHOCKS"
ptm <- proc.time()
uv.dists <- lapply(seq(ncol(u)), function(i) {
  student.t(nu = df, gamma = skew[i])
})
shocks <- dc.shocks(u, uv.dists, cluster)
proc.time() - ptm

"MARGINAL LIKELIHOODS"
ptm <- proc.time()
marginalLL <- marginalLogLikelihood(shocks, uv.dists, cluster)
marginalLL <- sum(marginalLL)
proc.time() - ptm

"JOINT LIKELIHOOD"
ptm <- proc.time()
jointLL <- jointLogLikelihood(
  shocks,
  uv.dists,
  mv.df = df,
  mv.skew = skew,
  alpha = alpha,
  beta = beta,
  cluster = cluster
)
proc.time() - ptm
jointLL

stopCluster(cluster)