rm(list = ls())

get_realized <- function(name) {
  load(sprintf('data/derived/crra/%s.RData', name))
  sapply(results, function(r) r$realized)
}

realized <- list(
  constant_norm = get_realized('simple_5_oos_constant_norm'),
  constant_std = get_realized('simple_5_oos_constant_std'),
  dynamic_norm = get_realized('simple_5_oos_dynamic_norm'),
  dynamic_std = get_realized('simple_5_oos_dynamic_std'),
  dynamic_ghst = get_realized('simple_5_oos_dynamic_ghst')
)

lapply(realized, mean)
lapply(realized, sd)
lapply(realized, function(r) mean(r) / sd(r))
lapply(realized, function(r) mean((1 + r) ^ (1 - 5)) ^ (1 / (1 - 5)) - 1)
