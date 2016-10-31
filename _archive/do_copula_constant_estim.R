#' ARCHIVED BECAUSE OLD COPULA INTERFACE AND GARCH
#' 
#' Estimate copulas with constant correlation matrices
#'
#' These estimations use the dynamic code, but with parameter restrictions
#' alpha = beta = 0.
#'

# THIS CODE NEEDS TO BE UPDATED TO WORK WITH NEW COPULA CODE!!!

# Workspace Setup --------------------------------------------------------
library(parallel)
library(devtools)
load_all('wimbledon')

rm(list = ls())

load('data/derived/garch_unires_model.RData')
u <- df.u[, -1]
rm(df.u)

optimizeFn <- function(theta, dist, data, cluster = NULL) {
  params <- dc.get.params(theta, ncol(data), dist)

  -wimbledon::dc.ll.total(
    dist.params = params$dist.params,
    alpha = params$alpha,
    beta = params$beta,
    u = data,
    cluster = cluster
  )
}

# Asymmetric Student's t Optimization ----
theta.ghskt <- c(
  11.80294,
  -0.05473225,
  0.07062225,
  -0.17046365,
  -0.12501372,
  0.09492990,
  0.02159156
)

ui <- rbind(
  c( 1, rep(0, ncol(u))), # min df
  c(-1, rep(0, ncol(u))), # max df

  cbind(0, diag( 1, ncol(u))), # min skew
  cbind(0, diag(-1, ncol(u)))  # max skew
)

ci <- cbind(c(
       6.0000,   # min df
     -20.0000,   # max df
  rep(-0.25, ncol(u)),  # min skew
  rep(-0.25, ncol(u))   # max skew
))

cluster <- prepare.cluster()

optim.ghskt <- constrOptim(
  theta.ghskt,
  optimizeFn,
  grad = NULL,

  dist = 'ghskt',
  data = u,
  cluster = cluster,

  ui = ui,
  ci = ci,

  control = list(
    trace = 6,
    maxit = 1000
  )
)

stopCluster(cluster)
model.copula.constant.ghskt <- build.output(u, optim.ghskt$par, 'ghskt')

save(
  model.copula.constant.ghskt,
  file = 'data/derived/model_copula_constant_ghskt.RData'
)

# Symmetric Student's t optimization ----
theta.ght <- c(11.80294)
cluster <- prepare.cluster()

optim.ght <- optim(
  theta.ght,
  optimizeFn,
  gr = NULL,

  dist = 'ght',
  data = u,
  cluster = cluster,

  method = "Brent",
  lower = 6,
  upper = 50
)

stopCluster(cluster)
model.copula.constant.ght <- build.output(u, optim.ght$par, 'ght')

save(
  model.copula.constant.ght,
  file = 'data/derived/model_copula_constant_ght.RData'
)

# Gaussian Copula ----
# No optimization to do here as correlation matrix is MoM estimated
model.copula.constant.gauss <- build.output(u, NULL, 'gauss')

save(
  model.copula.constant.gauss,
  file = 'data/derived/model_copula_constant_gauss.RData'
)
