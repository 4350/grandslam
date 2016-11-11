rm(list = ls())

library(dplyr)

MODEL_NAME <- 'results_full_dynamic_std_10000'
# MODEL_NAME <- 'results_full_sample'

STRATEGIES <- list(
  '5F',
  '5F_EXCL_HML',
  '5F_EXCL_CMA',
  '6F',
  '6F_EXCL_HML',
  '6F_EXCL_CMA'
)

load_field <- function(field) {
  out <- lapply(STRATEGIES, function(strategy) {
    load(sprintf('data/derived/mv/%s_%s.RData', MODEL_NAME, strategy))
    results[[field]]
  })
  
  names(out) <- STRATEGIES
  out
}

# Summary statistics of realized returns ----
returns <- load_field('portfolio_return')

stats <- lapply(returns, function(strategy) {
  stats_labels <- c('Mean', 'Stdev', 'Skewness', 'Kurtosis')
  stats <- as.list(fBasics::basicStats(strategy)[stats_labels, ])
  names(stats) <- stats_labels
  stats$Mean <- 52 * stats$Mean * 100
  stats$Stdev <- sqrt(52) * stats$Stdev * 100
  
  # SR
  SR <- stats$Mean / (stats$Stdev)
  MDD <- maxDrawdown(strategy)
  data.frame(stats, MDD = MDD, SR = SR)
})

cat(sprintf('== %s ==\n\n', MODEL_NAME))
print(bind_rows(stats, .id = 'Strategy'))

# Summary statistics of weights ----
weights <- load_field('weights')

wdf <- bind_rows(lapply(weights, function(weight) {
  # lol
  100 * data.frame(t(data.frame(apply(weight, 2, mean))))
}), .id = 'Strategy')

stargazer::stargazer(wdf, summary = FALSE, type = 'text', digits = 2)
