# XXX ----
library(ggplot2)
library(tidyr)
library(zoo)

get_portfolio_field <- function(field, files) {
  lapply(files, function(file) {
    load(sprintf('data/derived/mv/results_oos_%s.RData', file))
    results[[field]]
  })
}

weights <- get_portfolio_field('weights', list(
  'dynamic_ghst_10000_5F',
  'dynamic_ghst_10000_5F_EXCL_HML'
))

weights_df <- lapply(weights, function(model) {
  model <- apply(model, 2, function(w) rollmeanr(w, 26))
  
  model <- data.frame(
    t = 1:nrow(model),
    model
  )
  gather(model, 'Factor', 'Weights', -1)
})
names(weights_df) <- c('All', 'Excl. HML')

weights_df <- bind_rows(weights_df, .id = 'Strategy')
# weights_df <- lapply(weights,
#                      function(model) {
#                        data.frame(model) %>%
#                           gather(., 'Factor', 'Weights') %>%
#                           group_by(Factor) %>%
#                           rollmeanr(., 52)
#                      })



# 
# weights_df <- data.frame(
#   t = 1:(914 - 51),
#   HML_ALL = rollmeanr(weights[[1]][, 'HML'], 52),
#   HML = rollmeanr(weights[[2]][, 'HML'], 52)
# )

# gathered <- gather(weights_df, 'Strategy', 'Weights', -1)

ggplot(weights_df, aes(x = t, y = Weights, color = Strategy)) +
  facet_grid(. ~ Factor) +
  geom_line()
# ggplot(gathered, aes(x = t, y = Weights, color = Strategy)) +
#   geom_line()


# Returns ----------------------------------------------------------------

returns <- data.frame(
  t = 1:914,
  portfolio_returns(list(
    
  ))
)

returns_gathered <- gather(returns, Strategy, Return, -1)

ggplot(returns_gathered, aes(x = t, y = Return, color = Strategy)) +
  geom_line()

# Sharpe
apply(returns[, -1], 2, mean)
apply(returns[, -1], 2, sd)
apply(returns[, -1], 2, mean) / apply(returns[, -1], 2, sd) * 52 / sqrt(52)

