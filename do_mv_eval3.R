rm(list = ls())

library(dplyr)
library(PerformanceAnalytics)
library(stargazer)

MODEL_NAME <- 'results_full_dynamic_std_10000'
#MODEL_NAME <- 'results_full_sample'

STRATEGIES <- list(
  '5F',
  '5F_EXCL_CMA',
  '5F_EXCL_HML',
  '5F_EXCL_RMW',
  '6F',
  '6F_EXCL_CMA',
  '6F_EXCL_HML',
  '6F_EXCL_RMW'
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

stargazer::stargazer(t(wdf), summary = FALSE, type = 'text', digits = 3)

# T-testing diffs  ------------------------------------------------------------------------

t_test_factor <- function(strategy1, strategy2) {
  
  factors = c('Mkt.RF','SMB','Mom','HML','CMA','RMW')
  
  out.list <- lapply(factors, function(factor) {

    load(sprintf('data/derived/mv/%s_%s.RData', MODEL_NAME, strategy1))
    results1 <- results[['weights']]
    load(sprintf('data/derived/mv/%s_%s.RData', MODEL_NAME, strategy2)) 
    results2 <- results[['weights']]
    
    if(!(
      factor %in% colnames(results1) && factor %in% colnames(results2)
    )) {
      list()
    } else {
      
      weights1 <- results1[,factor] * 100
      weights2 <- results2[,factor] * 100
      
      df <- data.frame(weights1, weights2)
      colnames(df) <- c(strategy1, strategy2)
      
      test <- t.test(df[,strategy1], df[,strategy2], paired = T)
      
      list(
        factor = factor,
        strategy1 = strategy1,
        strategy2 = strategy2,
        statistic = round(test$statistic,2),
        estimate = round(test$estimate,1),
        se = round(test$estimate / test$statistic,4),
        p = round(test$p.value,4)
      )
      
    }
  })
  
  out.list <- bind_rows(out.list)
  
  stargazer(out.list, summary = F, type = 'text', digits = 2, digits.extra = 0)
  
}

t_test_factor('5F_EXCL_CMA', '5F')
t_test_factor('5F_EXCL_HML', '5F')


t_test_factor('6F_EXCL_CMA', '6F')
t_test_factor('6F_EXCL_HML', '6F')
