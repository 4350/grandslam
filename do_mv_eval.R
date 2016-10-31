library(ggplot2)
library(tidyr)

get_portfolio_field <- function(field, files) {
  lapply(files, function(file) {
    load(sprintf('data/derived/mv/results_oos_%s.RData', file))
    results[[field]]
  })
}

weights <- get_portfolio_field('weights', list(
  'constant_norm_100000_All',
  'indep_100000_All'
))

weights <- data.frame(
  t = 1:914,
  sapply(weights, function(w) w[, 1])
)

gathered <- gather(weights, Strategy, Market, -1)

ggplot(gathered, aes(x = t, y = Market, color = Strategy)) +
  geom_line()

# Returns ----------------------------------------------------------------

returns <- data.frame(
  t = 1:914,
  portfolio_returns(list(
    'dynamic_std_10000_All',
    'dynamic_norm_10000_All',
    'constant_std_10000_All',
    'constant_ghst_10000_All',
    'constant_norm_10000_All',
    'indep_100000_All'
  ))
)

returns_gathered <- gather(returns, Strategy, Return, -1)

ggplot(returns_gathered, aes(x = t, y = Return, color = Strategy)) +
  geom_line()

# Sharpe
apply(returns[, -1], 2, mean)
apply(returns[, -1], 2, sd)
apply(returns[, -1], 2, mean) / apply(returns[, -1], 2, sd) * 52 / sqrt(52)

