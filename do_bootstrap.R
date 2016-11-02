#' Run `do_bootstrap_index` first!

# Setup ------------------------------------------------------------------

rm(list = ls())

library(devtools)
library(doParallel)
library(tictoc)
load_all('australian')

load('data/derived/weekly-full.RData')
load('data/derived/garch/model_GARCH_chosen.RData')
load('data/derived/copula/full_dynamic.RData')
load('data/derived/bootstrap/bs_index.RData')

# Select relevant models. NULL out Omega in the fitted copula
MODEL_NAME <- 'ghst'
MODEL_COPULA <- dynamic_copula_fit$ghst
MODEL_COPULA_FIT_ARGS <- list(distribution = MODEL_NAME, constant = F)

# Change this to start bootstrap at a new index (for massively parallel
# operations)
START_INDEX <- 1

# Configure parallel
cl <- makePSOCKcluster(7)
clusterEvalQ(cl, library(devtools))
clusterEvalQ(cl, load_all('australian'))
registerDoParallel(cl)

# Use getspec so we DON'T set fixed pars, just extract the spec
MODEL_GARCH <- lapply(model.GARCH, getspec)

# NULL out Omega for fitting!
MODEL_COPULA$fit@dynamics@Omega <- NULL

# Do Bootstrap -----------------------------------------------------------

# All this code should be moved to a separate file

#' Fit best GARCH models to original data, returning ugarchfit objects
fit_garch <- function(garch_models, data) {
  N <- length(garch_models)
  
  lapply(seq(N), function(i) {
    ugarchfit(
      garch_models[[i]],
      data[, i],
      solver = "hybrid"
    )
  })
}

#' Extract uniforms from GARCH fits
get_uniforms <- function(garch_fits) {
  sapply(garch_fits, function(fit) {
    rugarch:::psghst(
      fit@fit$z,
      shape = fit@fit$coef[['shape']],
      skew = fit@fit$coef[['skew']]
    )
  })
}

dir.create(sprintf('data/derived/bootstrap/%s', MODEL_NAME), showWarnings = F)

for (b in seq(START_INDEX, ncol(BS_INDEX))) {
  tic(sprintf('Model "%s". Bootstrap %d', MODEL_NAME, b))
  
  # Need to make a data frame; something about tibbles messes up
  data <- data.frame(df[BS_INDEX[, b], -1])
  
  tic('GARCH fitting')
    fits <- fit_garch(MODEL_GARCH, data)
  toc()
  
  tic('Extract uniform residuals')
    uniforms <- get_uniforms(fits)
  toc()
  
  tic('Fit Copula')
  fitted_copula <- do.call(copula_fit, c(
    list(
      spec = MODEL_COPULA$fit,
      u = uniforms
    ),
    MODEL_COPULA_FIT_ARGS
  ))
  toc()
  
  results <- list(
    garch = fits,
    copula = fitted_copula$fit
  )
  
  save(
    results,
    file = sprintf('data/derived/bootstrap/%s/%d.RData', MODEL_NAME, b)
  )
  
  toc()
}