#' Bootstrap standard errors for our models
#'
#' This is done sequentially for each bootstrap, however, it is suitable
#' to run this process on multiple computers and later join the results.
#'
#' The script is written to output its results progressively, so you can cancel
#' it and the bootstrapped parameters so far will be available in an output
#' file.

# BOOTSTRAP PARAMETERS ---------------------------------------------------

library(devtools)
load_all('wimbledon')

rm(list = ls())

# RNG seed used for Bootstrap. Change this to run a parallel universe bootstrap
kRandomSeed <- 403

# Number of bootstrap repetitions to perform. Set this to a large number
kBSRepetitions <- 1000

# Bootstrap iteration to start at (<=kBootstrapRepetitions)
kBSIteration <- 4

# Average block length for stationary bootstrap
# See optim_block_length for choosing this number
kBSBlockLength <- 45

# Best models, as determined before
kGARCHModels <- list(
  Mkt.RF = garch.specgen(0, 0),
  HML = garch.specgen(1, 1),
  SMB = garch.specgen(1, 1),
  Mom = garch.specgen(1, 0),
  RMW = garch.specgen(1, 0),
  CMA = garch.specgen(1, 0)
)

# Output path
kBSOutputPath <- 'data/derived/bootstrap/ghskt'

# From previous runs, values around here are found to optimize the LL
load('data/derived/model_copula_dynamic_ghskt.RData')
kCopulaParams <- c(
  model.copula.dynamic.ghskt$dist.params$df,
  model.copula.dynamic.ghskt$dist.params$skew,
  model.copula.dynamic.ghskt$alpha,
  model.copula.dynamic.ghskt$beta
)

# c(
#   11.804937402,
#   -0.019263030,
#   0.064795233,
#   -0.161745377,
#   -0.148752424,
#   0.081236747,
#   0.004100415,
#   0.068927164,
#   0.912195831
# )

cz <- rep(0, 6)
kCopulaUi <- rbind(
  # Degrees of freedom bounds
  c( 1, cz, 0, 0),
  c(-1, cz, 0, 0),

  # Skewness Bounds
  cbind(cz, diag( 1, 6), cz, cz),
  cbind(cz, diag(-1, 6), cz, cz),

  # Alpha/Beta Bounds
  c(0, cz,  1,  0),
  c(0, cz,  0,  1),
  c(0, cz, -1, -1)
)
rm(cz)

kCopulaCi <- cbind(
  c(
      6.0000,       # min df
    -20.0000,       # max df (negative)
    rep(-0.25, 6),  # min skew
    rep(-0.25, 6),  # max skew (negative)
      0.0000,       # min alpha
      0.0000,       # min beta
     -0.9999        # max alpha + beta (negative)
  )
)

# Functions Unique to Each Copula Model ----

#' Function to optimize for GHSKT copula
#'
#' Unique because parameters are special for this copula
optimize <- function(params, u, cluster) {
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

#' Estimate GHSKT copula
#'
#' Unique because constraints need to be set specifically, and correct
#' optimized function must be called
estimate.copula <- function(u) {
  cluster <- makeCluster(detectCores() - 1)
  clusterEvalQ(cluster, library(ghyp))
  clusterEvalQ(cluster, library(devtools))
  clusterEvalQ(cluster, load_all('wimbledon'))

  param <- constrOptim(
    theta = kCopulaParams,
    f = optimize,
    grad = NULL,
    u = u,
    cluster = cluster,

    # Optimization constraints and control
    ui = kCopulaUi,
    ci = kCopulaCi,

    control = list(
      trace = 0,
      maxit = 1000
    )
  )

  stopCluster(cluster)
  param
}

#' Build bootstrap output for GHSKT copula
#'
#' Unique because names of parameters are hardcoded
build.bootstrap.output.copula <- function(b.copula) {
  b.out.copula <- c(
    b.copula$par,
    -b.copula$value,
    b.copula$convergence
  )
  names(b.out.copula) <- paste0('copula.', c(
    # Parameters
    'par.df',
    paste0('par.skew', seq(1, 6)),
    'par.alpha',
    'par.beta',

    # Diagnostics
    'll',
    'convergence'
  ))

  b.out.copula
}


# Bootstrap --------------------------------------------------------------
source('source/dynamic_bootstrap.R')

