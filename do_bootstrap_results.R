rm(list = ls())

library(australian)
library(dplyr)
library(stargazer)

# XXX Bootstrap should be named dynamic_norm or something
load('data/derived/copula/full_dynamic.RData')
PATH <- 'data/derived/bootstrap/norm'
estimate <- dynamic_copula_fit$norm

results <- lapply(list.files(PATH), function(data) {
  load(file.path(PATH, data))
  results
})

COPULA_PARAMETERS <- list(
  nu = function(c) c@distribution@nu,
  gamma1 = function(c) c@distribution@gamma[1],
  gamma2 = function(c) c@distribution@gamma[2],
  gamma3 = function(c) c@distribution@gamma[3],
  gamma4 = function(c) c@distribution@gamma[4],
  gamma5 = function(c) c@distribution@gamma[5],
  gamma6 = function(c) c@distribution@gamma[6],
  alpha = function(c) c@dynamics@alpha,
  beta = function(c) c@dynamics@beta
)

fit_results <- bind_rows(
  lapply(COPULA_PARAMETERS, function(param_fn) {
    param_se <- sd(sapply(results, function(r) param_fn(r$copula)))
    param_estimate <- param_fn(estimate$fit)
    names(param_estimate) <- NULL
    
    param_estimate[param_se == 0] <- NaN
    param_se[param_se == 0] <- NaN
    
    data.frame(
      estimate = param_estimate,
      se = param_se
    )
  }),
  .id = 'param'
)

fit_results$t <- fit_results$estimate / fit_results$se
fit_results$p <- 2 * (pt(-abs(fit_results$t), length(results)))

stargazer(
  fit_results,
  type = 'text',
  summary = F,
  digits = 4,
  digits.extra = 0
)
cat(sprintf('N: %d', length(results)))
