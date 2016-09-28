#' Estimate copula with constant correlation matrices

# Reset ----
rm(list = ls())
load('data/derived/garch_unires_model.RData')
u <- df.u[, -1]
rm(df.u)

library(ghyp)
library(parallel)
library(devtools)
load_all('wimbledon')

# Estimate Gaussian ----
load_all('wimbledon')

optimize.gauss <- function(params, u, cluster) {
  dist <- ghyp::gauss(
    mu = rep(0, ncol(u)),
    sigma = as.correlation.matrix(params)
  )
  
  wimbledon::cc.ll.total(u, dist, cluster)
}

#' Log likelihood of params with a Student t distribution
#'
#' @param params correlations, df, skewness
#' @param u uniform residuals
optimize.ghskt <- function(params, u, cluster) {
  num.corr <- ncol(u) * (ncol(u) - 1) / 2
  num.skew <- ncol(u)
  
  sigma <- as.correlation.matrix(params[1:num.corr])
  df <- params[num.corr + 1]
  
  skew <- rep(0, num.skew)
  if (length(params) > (num.corr + 1)) {
    skew <- params[(num.corr + 2):(num.corr + 1 + num.skew)]
  }
  
  dist <- ghyp::student.t(
    nu = df,
    mu = rep(0, ncol(u)),
    sigma = sigma,
    gamma = skew
  )
  
  wimbledon::cc.ll.total(u, dist, cluster)
}

prepare.cluster <- function() {
  cluster <- makeCluster(detectCores() - 1)
  clusterEvalQ(cluster, library(ghyp))
  clusterEvalQ(cluster, library(devtools))
  clusterEvalQ(cluster, load_all('wimbledon'))
  cluster
}

# optimize.gauss(rep(0.75, 15), u)
# optimize.ghskt(c(rep(0.75, 15), 7), u)

cluster <- prepare.cluster()
optimize.ghskt(
  c(rep(0.75, 15), 7, rep(0.90, ncol(u))),
  u = u,
  cluster = cluster
)

stopCluster(cluster)
rm(ptc)

# Estimate Gaussian ----
cluster <- prepare.cluster()
c <- cor(apply(u, 2, qnorm))
params <- c[lower.tri(c)]

param.gauss <- optim(
  params,
  optimize.gauss,
  u = u,
  cluster = cluster,
  control = list(
    trace = 6,
    fnscale = -1,
    maxit = 20000
  ),
  hessian = T
)

stopCluster(cluster)
rm(cluster)
