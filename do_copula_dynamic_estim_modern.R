# Reset ----

rm(list = ls())
load('data/derived/garch_unires_model.RData')
u <- df.u[, -1]
rm(df.u)

# Select factors
u <- u[,c('Mkt.RF','SMB','Mom','RMW','CMA')]

library(ghyp)
library(parallel)
library(devtools)
load_all('wimbledon')

prepare.cluster <- function() {
  cluster <- makeCluster(detectCores() - 1)
  clusterEvalQ(cluster, library(ghyp))
  clusterEvalQ(cluster, library(devtools))
  clusterEvalQ(cluster, load_all('wimbledon'))
  cluster
}

# Estimate Asymmetric T-distribution ----
optimize.ghskt <- function(params, u, cluster) {
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

cluster <- prepare.cluster()

# From previous run
params <- c(
  11.804937402,
  -0.019263030,
  -0.161745377,
  -0.148752424,
   0.081236747,
   0.004100415,
   0.068927164,
   0.912195831
)

param.ghskt <- constrOptim(
  theta = params,
  optimize.ghskt,
  grad = NULL,
  u = u,
  cluster = cluster,
  ui = rbind(
    c( 1,  0,  0,  0,  0,  0,  0,  0),
    c(-1,  0,  0,  0,  0,  0,  0,  0),
    
    c( 0,  1,  0,  0,  0,  0,  0,  0),
    c( 0, -1,  0,  0,  0,  0,  0,  0),
    c( 0,  0,  1,  0,  0,  0,  0,  0),
    c( 0,  0, -1,  0,  0,  0,  0,  0),
    c( 0,  0,  0,  1,  0,  0,  0,  0),
    c( 0,  0,  0, -1,  0,  0,  0,  0),
    c( 0,  0,  0,  0,  1,  0,  0,  0),
    c( 0,  0,  0,  0, -1,  0,  0,  0),
    c( 0,  0,  0,  0,  0,  1,  0,  0),
    c( 0,  0,  0,  0,  0, -1,  0,  0),
    
    c( 0,  0,  0,  0,  0,  0,  1,  0),
    c( 0,  0,  0,  0,  0,  0,  0,  1),
    c( 0,  0,  0,  0,  0,  0, -1, -1)
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
rm(cluster)
rm(optimize.ghskt)
rm(params)
save(
  param.ghskt,
  file = ('data/derived/dcopula_param_ghskt_modern.RData')
)