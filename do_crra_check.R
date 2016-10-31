rm(list = ls())

get_realized <- function(name) {
  load(sprintf('data/derived/crra/%s.RData', name))
  sapply(results, function(r) r$realized)
}

realized <- list(
  constant_std = get_realized('simple_5_full_constant_std'),
  dynamic_ghst = get_realized('simple_5_full_dynamic_ghst'),
  indep = get_realized('simple_5_full_indep')
)

lapply(realized, mean)
lapply(realized, sd)
lapply(realized, function(r) mean(r) / sd(r))
lapply(realized, function(r) mean((1 + r) ^ (1 - 5)) ^ (1 / (1 - 5)) - 1)
