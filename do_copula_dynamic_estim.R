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

estimate.dynamic.copula <- function(theta, dist) {
  # Closes over "data"
  data <- u

  cluster <- prepare.cluster()

  # Constraints for the optimizers
  constr <- dc.constraints(ncol(u))
  ui <- constr$ui$alphabeta
  ci <- constr$ci$alphabeta
  if (dist != 'gauss') {
    if (dist == 'ghskt') {
      ui <- adiag(constr$ui$df, constr$ui$skew, ui)
      ci <- c(constr$ci$df, constr$ci$skew, ci)
    } else {
      ui <- adiag(constr$ui$df, ui)
      ci <- c(constr$ci$df, ci)
    }
  }

  result <- constrOptim(
    theta,
    dc.optimize.fn,
    grad = NULL,

    dist = dist,
    data = data,
    cluster = cluster,

    ui = ui,
    ci = ci,

    control = list(
      trace = 6,
      maxit = 1000
    )
  )

  stopCluster(cluster)

  # Respond with full model thing
  build.output(u, result$par, dist)
}

# Estimate Asymmetric T-distribution ----
model.copula.dynamic.ghskt <- estimate.dynamic.copula(
  theta = c(
    11.804937402,
    -0.019263030,
     0.064795233,
    -0.161745377,
    -0.148752424,
     0.081236747,
     0.004100415,
     0.068927164,
     0.912195831
  ),
  dist = 'ghskt'
)
save(
  model.copula.dynamic.ghskt,
  file = 'data/derived/model_copula_dynamic_ghskt.RData'
)

# Estimate Symmetric T-distribution ----
model.copula.dynamic.ght <- estimate.dynamic.copula(
  theta = c(
    11.78368353,
     0.06881473,
     0.91247336
  ),
  dist = 'ght'
)
save(
  model.copula.dynamic.ght,
  file = 'data/derived/model_copula_dynamic_ght.RData'
)

# Estimate Gaussian Copula -----------------------------------------------
model.copula.dynamic.gauss <- estimate.dynamic.copula(
  theta = c(
    0.06577153,
    0.9145452
  ),
  dist = 'gauss'
)
save(
  model.copula.dynamic.gauss,
  file = 'data/derived/model_copula_dynamic_gauss.RData'
)

# Information criteria for parameters ----

# CODE BELOW DOES NOT WORK WITH NEW COPULA MODEL STRUCTURE
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

