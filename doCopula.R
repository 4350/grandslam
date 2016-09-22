# Estimate Copula parameters on GARCH residuals

# Libraries ----
library(ghyp)
library(dplyr)
library(parallel)
library(compiler)


# Reset Workspace ----
rm(list = ls())
load('data/derived/GARCHunifres.RData')
u <- dfWeekly.u[, -1]
rm(dfWeekly.u)
# load('data/derived/GARCHunifres.RData')
# u <- as.matrix(dfDaily.u[, -1])
# rm(dfDaily.u)
# load('data/derived/garch5f-u.Rdata')

# Optimization Boilerplate ----

enableJIT(3)
cmpfile('func/dynamicCopula.R')
cmpfile('func/dghsstx.R')
loadcmp('func/dghsstx.Rc')
loadcmp('func/dynamicCopula.Rc')


prepareCluster <- function() {
  cluster <- makeCluster(detectCores() - 1)
  clusterEvalQ(cluster, library(ghyp))
  clusterExport(cluster, "marginalLogLikelihood")
  clusterExport(cluster, "jointLogLikelihood")
  clusterExport(cluster, "dc.shocks")
  clusterExport(cluster, "dc.shocks.std")
  clusterExport(cluster, "dc.Q")
  clusterExport(cluster, "dc.Omega")
  clusterExport(cluster, "dc.Correlation")
  clusterExport(cluster, "besselM3")
  clusterExport(cluster, "dghsstx")
  cluster
}

#' Optimize log-likelihood with free skewness parameters
optimize2Fn <- function(params, u, cluster) {
  df <- params[1]
  skew <- params[2:(ncol(u) + 1)]
  alpha <- tail(params, 2)[1]
  beta <- tail(params, 1)[1]
  
  -totalLogLikelihood(u, df, skew, alpha, beta, cluster)
}

optimize1Fn <- function(params, u, cluster) {
  df <- params[1]
  skew <- rep(params[2], ncol(u))
  alpha <- params[3]
  beta <- params[4]
  
  print(params)
  
  -totalLogLikelihood(u, df, skew, alpha, beta, cluster)
}

# Optimization with single Lambda ----
params <- c(10.74819641, -0.02613533,  0.06781287,  0.91311380)
#params <- c(9, 0, 0.03, 0.96)
cluster <- prepareCluster()

ptc <- proc.time()
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
      8.00,
    -12.00,
     -0.15,
     -0.15,
      0.00,
      0.00,
     -0.9999
  ),
  control = list(
    trace = 6,
    maxit = 1000
  ),
  method = 'Nelder-Mead'
)

proc.time() - ptc

stopCluster(cluster)


# Optimization Runner ----

# These parameters break qghyp if one uses a low tolerance because there are
# no roots in the intervals it looks at...
#
# Note: The specialized multivariate density function does not work for
# symmetric t; so don't try gamma = 0
params <- c(15, rep(0.10, 6), 0.03, 0.96)

# Best params as of 2016-09-21 for df (8-12), gamma (-0.20-0.20)
params <- c(
  10.796258944,
  0.035832769,
  -0.184606557,
  0.009726851,
  0.084787508,
  -0.036544851,
  -0.152462968,
  0.067887334,
  0.912931257
)
# params <- c(5.53417002054087, 0.0218613978265847, -0.0141872411754378, 
#             -0.133528167821921, -0.0160852945932761, 0.0621979871082812, 
#             0.00687898447016727, 0.0378971855489661, 0.956425966391781)
#params <- c(5, rep(0, ncol(u)), 0.03, 0.96)
cluster <- prepareCluster()

ptc <- proc.time()
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
    6,       # min df
    -20,     # -(max df)
    
    # SMB daily cannot handle skew = 0.50 -- which seems like a large value
    # from looking at Christoffersen; constrain hard
    -0.25,      # min skew
    -0.25,      # -(max skew)
    -0.25,      # min skew
    -0.25,      # -(max skew)
    -0.25,      # min skew
    -0.25,      # -(max skew)
    -0.25,      # min skew
    -0.25,      # -(max skew)
    -0.25,      # min skew
    -0.25,      # -(max skew)
    -0.25,      # min skew
    -0.25,      # -(max skew)
    
    0,       # min alpha
    0,       # min beta
    -0.9999  # -(max alpha + beta)
  ),
  control = list(
    trace = 6,
    maxit = 1000
  )
)
proc.time() - ptc

stopCluster(cluster)

# Diagnostics ----
