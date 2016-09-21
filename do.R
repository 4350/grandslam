#' Estimate a GARCH-COPULA model on Weekly Returns Data a la
#' Christoffersen and Langlois (2013)

# Libraries ----
library(rugarch)
library(sn)
library(ghyp)
library(dplyr)
library(parallel)
library(compiler)

# Reset Workspace to just Returns Data ----
rm(list = ls())
load("data/ff-weekly.RData")
df <- select(df, -RF)

# Build spec and fit GARCH model ----
garch.spec <- ugarchspec(
  mean.model = list(
    armaOrder = c(3, 0)
  ),
  variance.model = list(
    model = "fGARCH",
    submodel = "NGARCH",
    variance.targeting = TRUE
  ),
  distribution.model = 'sstd'
)

# Fit GARCH to each factor, storing result in a list for each factor,
# where each factor in turn has a list with element "fit" holding the results
# (so we can hold more stuff next to it later)
garch.fit <- apply(
  select(df, -Date), 2,
  function(ts) {
    list(fit = ugarchfit(garch.spec, ts))
  }
)

# Find the eta elements for the copula from the standardized residuals ----
# d = pdf, p = cdf, q = quantile

get.u <- function(ts) {
  residuals.st <- ts$fit@fit$residuals / ts$fit@fit$sigma
  
  # Generate empirical CDF for the standardized residuals and use it to
  # compute union margins (eta)
  #cdf <- ecdf(residuals.st)
  
  # Ugly hack. PS get your shit together Victor
  cdf <- rugarch:::psstd(
    residuals.st,
    skew = ts$fit@fit$coef[['skew']],
    shape = ts$fit@fit$coef[['shape']]
  )
  
  list(
    fit = ts$fit,
    residuals.st = residuals.st,
    u = cdf#(residuals.st)
  )
}

garch.fit <- lapply(garch.fit, get.u)
u <- sapply(garch.fit, function(fit) fit$u)
save(u, file = 'data/derived/garch5f-u.Rdata')

# Dynamic Copula Estimation ----

enableJIT(3)
cmpfile('func/dynamicCopula.R')
loadcmp('func/dynamicCopula.Rc')

#' Optimize log-likelihood with a single skewness parameter
#'
#' @param params (df, skew, alpha, beta)
#' @param u 
#'
#' @return negative log likelihood
optimize1Fn <- function(params, u, cluster) {
  df <- params[1]
  skew <- rep(params[2], ncol(u))
  alpha <- params[3]
  beta <- params[4]
  
  -totalLogLikelihood(u, df, skew, alpha, beta, cluster)
}

#' Optimize log-likelihood with free skewness parameters
#'
#' @param params
#' @param u
#' @param cluster
#'
#' @return
#' @export
#'
#' @examples
optimize2Fn <- function(params, u, cluster) {
  df <- params[1]
  skew <- params[2:(ncol(u) + 1)]
  alpha <- tail(params, 2)[1]
  beta <- tail(params, 1)[1]
  
  -totalLogLikelihood(u, df, skew, alpha, beta, cluster)
}

#' Optimize conditional log-likelihood with a single skewness parameter
optimize1FnCond <- function(params, u, cluster) {
  df <- params[1]
  skew <- rep(params[2], ncol(u))
  alpha <- params[3]
  beta <- params[4]
  
  -condLogLikelihood(u, df, skew, alpha, beta, cluster)
}

# TESTING OF LOG LIKELIHOOD FNs ----
# params <- c(5, 1, 0.03, 0.96)
# u <- sapply(garch.fit, function(f) f$u)
# 
# kNumCores <- detectCores() - 1
# cluster <- makeCluster(kNumCores)
# clusterEvalQ(cluster, library(ghyp))
# clusterExport(cluster, "marginalLogLikelihood")
# clusterExport(cluster, "jointLogLikelihood")
# clusterExport(cluster, "dc.shocks")
# clusterExport(cluster, "dc.shocks.std")
# clusterExport(cluster, "dc.Q")
# clusterExport(cluster, "dc.Omega")
# clusterExport(cluster, "dc.Correlation")
# 
# # Start the clock!
# ptm <- proc.time()
# optimize1FnCond(params, u, cluster)
# proc.time() - ptm
# 
# ptm <- proc.time()
# optimize1Fn(params, u, cluster)
# proc.time() - ptm
# 
# stopCluster(cluster)

# Do Optimization ----
# Input data
params <- c(5, 0.05, 0.03, 0.96)
load('data/derived/garch5f-u.RData')

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

param.constrOptim <- constrOptim(
  theta = params,
  optimize1Fn,
  grad = NULL,
  u = u,
  cluster = cluster,
  ui = rbind(
    c( 1,  0,  0,  0),
    c(-1,  0,  0,  0),
    c( 0,  1,  0,  0),
    c( 0, -1,  0,  0),
    c( 0,  0,  1,  0),
    c( 0,  0,  0,  1),
    c( 0,  0, -1, -1)
  ),
  ci = rbind(
      4,     # min df
    -10,     # -(max df)
     -1,     # min skew
     -1,     # -(max skew)
      0,     # min alpha
      0,     # min beta
     -0.9999 # -(max alpha + beta)
  ),
  control = list(
    trace = 5,
    maxit = 100
  )
)

stopCluster(cluster)

# Perform optimization with free skewness per N ----

load('data/derived/garch5f-u.RData')
params <- c(5, rep(0, ncol(u)), 0.03, 0.96)

# Prepare clusters
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

param.constrOptim <- constrOptim(
  theta = params,
  optimize2Fn,
  grad = NULL,
  u = u,
  cluster = cluster,
  ui = rbind(
    c( 1,  0,  0,  0,  0,  0,  0,  0,  0),
    c(-1,  0,  0,  0,  0,  0,  0,  0,  0),
    
    c( 0,  1,  0,  0,  0,  0,  0,  0,  0),
    c( 0, -1,  0,  0,  0,  0,  0,  0,  0),
    c( 0,  0,  1,  0,  0,  0,  0,  0,  0),
    c( 0,  0, -1,  0,  0,  0,  0,  0,  0),
    c( 0,  0,  0,  1,  0,  0,  0,  0,  0),
    c( 0,  0,  0, -1,  0,  0,  0,  0,  0),
    c( 0,  0,  0,  0,  1,  0,  0,  0,  0),
    c( 0,  0,  0,  0, -1,  0,  0,  0,  0),
    c( 0,  0,  0,  0,  0,  1,  0,  0,  0),
    c( 0,  0,  0,  0,  0, -1,  0,  0,  0),
    c( 0,  0,  0,  0,  0,  0,  1,  0,  0),
    c( 0,  0,  0,  0,  0,  0, -1,  0,  0),
    
    c( 0,  0,  0,  0,  0,  0,  0,  1,  0),
    c( 0,  0,  0,  0,  0,  0,  0,  0,  1),
    c( 0,  0,  0,  0,  0,  0,  0, -1, -1)
  ),
  ci = rbind(
    4,       # min df
    -10,     # -(max df)
    
    -1,      # min skew
    -1,      # -(max skew)
    -1,      # min skew
    -1,      # -(max skew)
    -1,      # min skew
    -1,      # -(max skew)
    -1,      # min skew
    -1,      # -(max skew)
    -1,      # min skew
    -1,      # -(max skew)
    -1,      # min skew
    -1,      # -(max skew)
    
    0,       # min alpha
    0,       # min beta
    -0.9999  # -(max alpha + beta)
  ),
  control = list(
    trace = 6
  )
)

stopCluster(cluster)