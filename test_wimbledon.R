# Reset ----

rm(list = ls())
load('data/derived/GARCHunifres.RData')
u <- dfWeekly.u[, -1]
rm(dfWeekly.u)

library(ghyp)
library(parallel)
load_all('wimbledon')

prepare.cluster <- function() {
  cluster <- makeCluster(detectCores() - 1)
  clusterEvalQ(cluster, library(ghyp))
  clusterEvalQ(cluster, library(devtools))
  clusterEvalQ(cluster, load_all('wimbledon'))
  cluster
}

optimize.2 <- function(params, u, cluster) {
  df <- params[1]
  skew <- params[2:(ncol(u) + 1)]
  alpha <- tail(params, 2)[1]
  beta <- tail(params, 1)[1]
  
  -wimbledon::dc.ll.total(
    u,
    list(
      df = df,
      skew = skew
    ),
    alpha = alpha,
    beta = beta,
    cluster
  )
}

# Estimate ----
cluster <- prepare.cluster()

params <- c(
  11.791643017,
  -0.018823612,
   0.064814302,
  -0.161740839,
  -0.148497365,
   0.004072006,
   0.080768733,
   0.068956404
)

param.constrOptim <- constrOptim(
  theta = params,
  optimize.2,
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

stopCluster(cluster)
