#do_cdb_mv_tables


# Setup -------------------------------------------------------------------

rm(list = ls())

library(dplyr)
library(PerformanceAnalytics)
library(stargazer)

#MV
#MODEL_NAME <- 'results_full_dynamic_std_10000'
#MODEL_NAME <- 'results_sample'

#CDB
MODEL_NAME <- 'full_dynamic_std_10000'

STRATEGIES <- list(
  '5F',
  '5F_EXCL_HML',
  '5F_EXCL_CMA',
  '5F_EXCL_RMW',
  '6F',
  '6F_EXCL_HML',
  '6F_EXCL_CMA',
  '6F_EXCL_RMW'
)

load_field <- function(field) {
  out <- lapply(STRATEGIES, function(strategy) {
    if(MODEL_NAME == 'full_dynamic_std_10000') {
      load(sprintf('data/derived/cdb/constrOptim_q5_%s_%s.RData', MODEL_NAME, strategy))
    } else {
      load(sprintf('data/derived/mv/%s_%s.RData', MODEL_NAME, strategy))
    }
    
    results[[field]]
  })
  
  names(out) <- STRATEGIES
  out
}



# Summary stats  ----------------------------------------------------------

returns <- load_field('portfolio_return')
cdb <- load_field('cdb')
var <- load_field('var')

stats <- lapply(returns, function(strategy) {
  stats_labels <- c('Mean', 'Stdev', 'Skewness', 'Kurtosis')
  stats <- as.list(fBasics::basicStats(strategy)[stats_labels, ])
  names(stats) <- stats_labels
  stats$Mean <- 52 * stats$Mean * 100
  stats$Stdev <- sqrt(52) * stats$Stdev * 100
  
  # SR
  SR <- stats$Mean / (stats$Stdev)
  MDD <- maxDrawdown(strategy) * 100
  data.frame(stats, SR = SR, MDD = MDD)
})

var_avg <- lapply(var, function(strategy) mean(strategy))
cdb_avg <- lapply(cdb, function(strategy) mean(strategy))

out_stats <- t(
  cbind(
    bind_rows(stats),
    VaR = unlist(var_avg),
    CDB = unlist(cdb_avg)
  )
)


# Optimal weights ---------------------------------------------------------

weights <- load_field('weights')

wdf <- t(bind_rows(lapply(weights, function(weight) {
  # lol
  100 * data.frame(t(data.frame(apply(weight, 2, mean))))
})))

# Output ------------------------------------------------------------------

out_matrix <- rbind(
  wdf, 
  out_stats
)

colnames(out_matrix) <- STRATEGIES
stargazer(out_matrix, summary = FALSE, type = 'text', digits = 2)
#stargazer(out_matrix, summary = FALSE, type = 'latex', digits = 2, header = FALSE)

