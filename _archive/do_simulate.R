#' Simulate Return Series from ARMA-GARCH-DCC Copula Model
#'
#' For each T,
#'
#' Generate a battery of shocks from copula MV distribution using
#' conditional correlation matrix (for t = 0, use unconditional)
#'
#' Then, to update correlation matrix for next period:
#'
#' 1. Standardize shocks
#' 2. Send through DCC model to update correlation matrix for next period
#'
#' To generate return series:
#'
#' 1. Perform uniform transformation of (non-standardized) shock according to
#'    copula univariate distributions
#' 2. Compute quantile according to each GARCH model to get GARCH shocks
#' 3. Send through ARMA-GARCH model to get return for this period
#'
#' Nothing in the copula depends on the GARCH steps; it makes sense to generate
#' the full copula series of shocks first, and then transform these shocks into
#' returns.
#'
#' The copula simulation is hard to run in parallel, as it involves
#' multivariate recursion.
#'
#' In order to simulate, the models must be fully specified.
#'
#' For the copula, I need this:
#'
#'      list(
#'        dist.params = list(
#'          df = 10,             # NULL if normal
#'          skew = c(3, 5, 6, 3) # NULL if symmetric
#'        ),
#'        alpha = 0.06,
#'        beta = 0.91,
#'        Omega = matrix(c(
#'          1.00, 0.50, 0.50, 0.50,
#'          0.50, 1.00, 0.50, 0.50,
#'          0.50, 0.50, 1.00, 0.50,
#'          0.50, 0.50, 0.50, 1.00
#'        ), nrow = 4, byrow = T)
#'      )
#'
#'  For each ARMA-GARCH, I need the spec and parameters. This should work
#'  fine with whatever input/output rugarch deals with.

# Setup Libraries ----
library(devtools)
library(rugarch)
library(tictoc)
library(uuid)
library(parallel)
load_all('wimbledon')

# SETUP ------------------------------------------------------------------
rm(list = ls())

load('data/derived/model_GARCH.RData')
load('data/derived/model_copula_dynamic_ghskt.RData')

# These are the models that will be used for simulation. The format is
# ad hoc for the copula model, but matches output from estimation.
#
# The GARCH models should be a list of uGARCHfit objects (ugarchspec doesn't
# work because of reasons)
kCopulaModel <- model.copula.dynamic.ghskt
kGARCHModels <- model.GARCH
kOutputPath <- 'data/derived/simulation/ghskt'

# Number of simulations per run
T <- 4500

# Parallel simulations: Give number of simulations to run at the same time
B <- parallel::detectCores() - 1

rm(model.GARCH, model.copula.dynamic.ghskt)

# COPULA SIMULATION LOOP -------------------------------------------------

dir.create(kOutputPath, recursive = T, showWarnings = F)

cluster <- makeCluster(B)
clusterEvalQ(cluster, library(devtools))
clusterEvalQ(cluster, library(ghyp))
clusterEvalQ(cluster, library(rugarch))
clusterEvalQ(cluster, load_all('wimbledon'))

simulate <- function(i, T, garch, copula) {
  # Simulate
  result <- sim.c(T, garch, copula)

  # Pretty results
  data.frame(
    series = result$series,
    sigma = result$sigma,
    resid = result$resid
  )
}

n.simulations = 0
while (TRUE) {
  tic()
  results <- parLapplyLB(
    cluster,
    seq(B),
    simulate,
    T = T,
    garch = kGARCHModels,
    copula = kCopulaModel
  )
  toc()

  # Save to unique file
  lapply(results, function(result) {
    file <- file.path(kOutputPath, paste0(UUIDgenerate(), '.csv'))
    write.csv(result, file = file)
  })

  n.simulations <- n.simulations + length(results)
  cat('Total: ', n.simulations, '\n')
}

stopCluster(cluster)
