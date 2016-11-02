
# Setup ------------------------------------------------------------------

rm(list = ls())

library(dplyr)
library(ggplot2)

threshold_correlations_pairs <- function(stdresid, pairs) {
  results <- lapply(pairs, function(pair) {
    # qs really should be an argument to correlations.threshold
    coef <- correlations.threshold(stdresid, pair[1], pair[2])['coef', ]
    qs <- seq(0.10, 0.90, by=0.01)
    
    data.frame(
      factor1 = pair[1],
      factor2 = pair[2],
      qs = qs,
      coef = coef
    )
  })
  
  bind_rows(results)
}

threshold_correlations_models <- function(models, pairs) {
  results <- lapply(models, function(model) {
    load(sprintf('data/derived/stdresid/full_constant_%s.RData', model))
    threshold_correlations_pairs(stdresid, pairs)
  })
  
  names(results) <- models
  
  bind_rows(results, .id = 'model')
}

# Go go go ---------------------------------------------------------------

PAIRS <- list(
  c('Mkt.RF', 'HML'),
  c('Mkt.RF', 'RMW'),
  c('Mkt.RF', 'CMA'),
  c('Mom', 'HML'),
  c('Mom', 'RMW'),
  c('Mom', 'CMA'),
  c('CMA', 'HML'),
  c('RMW', 'HML'),
  c('RMW', 'CMA')
)

MODELS <- list(
  'norm',
  'std',
  'ghst'
)

tic(sprintf(
  'Threshold correlations for %d pairs, %d models',
  length(PAIRS),
  length(MODELS))
)
models <- threshold_correlations_models(MODELS, PAIRS)
toc()

save(models, file = 'data/derived/correlations_threshold_copula.RData')

# Plotting ---------------------------------------------------------------

ggplot(models, aes(x = qs, y = coef, color = model)) +
  geom_line() +
  facet_grid(factor1 ~ factor2)