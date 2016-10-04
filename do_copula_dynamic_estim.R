#' Estimate copulas with dynamic correlations

# Workspace Setup --------------------------------------------------------
library(parallel)
library(magic)
library(devtools)
load_all('wimbledon')

rm(list = ls())

load('data/derived/garch_unires_model.RData')
u <- df.u[, -1]
rm(df.u)

# Constraints for the optimizers
constr <- dc.constraints(ncol(u))

# Estimate Asymmetric T-distribution ----
cluster <- prepare.cluster()

# From previous run
theta.ghskt <- c(
  11.804937402,
  -0.019263030,
   0.064795233,
  -0.161745377,
  -0.148752424,
   0.081236747,
   0.004100415,
   0.068927164,
   0.912195831
)

optim.ghskt <- constrOptim(
  theta.ghskt,
  dc.optimize.fn,
  grad = NULL,

  dist = 'ghskt',
  data = u,
  cluster = cluster,

  ui = adiag(
    constr$ui$df,
    constr$ui$skew,
    constr$ui$alphabeta
  ),
  ci = c(
    constr$ci$df,
    constr$ci$skew,
    constr$ci$alphabeta
  ),

  control = list(
    trace = 6,
    maxit = 1000
  )
)

stopCluster(cluster)

model.copula.dynamic.ghskt <- build.output(
  u,
  dc.get.params(optim.ghskt$par, ncol(u), 'ghskt')
)
save(
  model.copula.dynamic.ghskt,
  file = ('data/derived/model_copula_dynamic_ghskt.RData')
)

# Estimate Symmetric T-distribution ----

# DO NOT RUN
#
# CODE BELOW DOES NOT SAVE ESTIMATES CORRECTLY
# SEE CODE ABOVE FOR CORRECT POST-FIDDLING TO SAVE EVERYTHING
#

optimize.ght <- function(params, u, cluster) {
  df <- params[1]
  alpha <- params[2]
  beta <- params[3]

  -wimbledon::dc.ll.total(
    u,
    list(
      df = df
    ),
    alpha = alpha,
    beta = beta,
    cluster
  )
}

params <- c(
  10,
  0.03,
  0.96
)

cluster <- prepare.cluster()

param.ght <- constrOptim(
  params,
  optimize.ght,
  grad = NULL,
  u = u,
  cluster = cluster,
  ui = rbind(
    c(  1,  0,  0),
    c( -1,  0,  0),
    c(  0,  1,  0),
    c(  0,  0,  1),
    c(  0, -1, -1)
  ),
  ci = rbind(
      4.00,  # min df (finite variance)
    -20.00,  # max df (arbitrary)
      0.00,  # min alpha
      0.00,  # min beta
     -0.9999 # max alpha + beta
  ),
  control = list(
    trace = 6,
    maxit = 1000
  )
)

stopCluster(cluster)
rm(optimize.ght)
rm(cluster)
rm(params)

model.copula.dynamic.ght <- param.ght
save(
  model.copula.dynamic.ght,
  file = ('data/derived/model_copula_dynamic_ght.RData')
)

# Estimate Gaussian Copula ----
optimize.gauss <- function(params, u, cluster) {
  alpha <- params[1]
  beta <- params[2]

  -wimbledon::dc.ll.total(
    u,
    list(),
    alpha = alpha,
    beta = beta,
    cluster
  )
}

params <- c(0.03, 0.96)
cluster <- prepare.cluster()

param.gauss <- constrOptim(
  params,
  optimize.gauss,
  grad = NULL,
  u = u,
  cluster = cluster,
  ui = rbind(
    c(  1,  0),
    c(  0,  1),
    c( -1, -1)
  ),
  ci = rbind(
    0.00,   # min alpha
    0.00,   # min beta
    -0.9999 # max alpha + beta
  ),
  control = list(
    trace = 6,
    maxit = 1000
  )
)

stopCluster(cluster)
rm(optimize.gauss)
rm(cluster)
rm(params)

model.copula.dynamic.gauss <- param.gauss
save(
  model.copula.dynamic.gauss,
  file = ('data/derived/model_copula_dynamic_gauss.RData')
)

# Information criteria for parameters ----
load('data/derived/model_copula_dynamic_ghskt.RData')
load('data/derived/model_copula_dynamic_ght.RData')
load('data/derived/model_copula_dynamic_gauss.RData')
load('data/derived/garch_unires_model.RData')
u <- df.u[, -1]
rm(df.u)

bic <- function(param) {
  -2 * -param$value + length(param$par) * log(nrow(u))
}
aic <- function(param) {
  -2 * -param$value + 2 * length(param$par)
}

# Increasing order of model complexity
param <- list(
  "Gaussian" = model.copula.dynamic.gauss,
  "Symmetric" = model.copula.dynamic.ght,
  "Skewed" = model.copula.dynamic.ghskt
)

models <- rbind(
  sapply(param, function(p) -p$value),
  sapply(param, bic),
  sapply(param, aic),

  # Christoffersen appears to count correlations estimated for Q
  sapply(param, function(p) length(p$par))
)
rownames(models) <- c('Log-l', 'BIC', 'AIC', 'Params')
colnames(models) <- names(param)
models

