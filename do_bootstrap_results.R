rm(list = ls())

library(australian)
library(dplyr)
library(stargazer)

copula_results <- function(estimate, bootstrap) {
  # Load bootstrap results from file
  path <- file.path('data/derived/bootstrap', bootstrap)
  results <- lapply(list.files(path), function(file) {
    load(file.path(path, file))
    results
  })
  
  # parameters we're interested in
  copula_parameters <- list(
    nu = function(c) c@distribution@nu,
    nu_rep = function(c) 1/c@distribution@nu,
    gamma1 = function(c) c@distribution@gamma[1],
    gamma2 = function(c) c@distribution@gamma[2],
    gamma3 = function(c) c@distribution@gamma[3],
    gamma4 = function(c) c@distribution@gamma[4],
    gamma5 = function(c) c@distribution@gamma[5],
    gamma6 = function(c) c@distribution@gamma[6],
    alpha = function(c) c@dynamics@alpha,
    beta = function(c) c@dynamics@beta
  )
  
  params <- lapply(copula_parameters, function(fn) {
    # Extract the parameter estimate from the fitted copula
    coef <- fn(estimate$fit)
    names(coef) <- NULL
    
    # For each bootstrap result, extract the copula parameter
    # and then compute the standard deviation of those parameter estimates
    se <- sd(sapply(results, function(r) fn(r$copula)))
    
    # Don't include unused parameters in nested models
    coef[se == 0 || is.nan(se)] <- NaN
    se[se == 0] <- NaN
    
    t <- coef / se
    p <- 2 * (pt(-(abs(t)), length(results)))
    data.frame(coef, se, p)
  })

  
  list(
    fit = bind_rows(params, .id = 'param'),
    N = length(results)
  )
}

load('data/derived/copula/full_dynamic.RData')
load('data/derived/copula/full_constant.RData')

stargazer(
  copula_results(constant_copula_fit$ghst, 'ghst')$fit,
  type = 'text',
  digits = 2,
  digits.extra = 0,
  summary = F,
  rownames = F
)

# lapply(dynamic_copula_fit, function(r) r$fit@dynamics@alpha + r$fit@dynamics@beta)
lapply(dynamic_copula_fit, function(r) r$ll)
