library(ggplot)
library(tidyr)

rm(list = ls())
DISTRIBUTION_NAME = 'ghskt'

load(sprintf('data/derived/var_es_constant_%s.RData', DISTRIBUTION_NAME))
var_es_constant <- var_es

load(sprintf('data/derived/var_es_dynamic_%s.RData', DISTRIBUTION_NAME))
var_es_dynamic <- var_es

load('data/derived/weekly-full.RData')

# Cut of df to VaR sample period
df <- df[(nrow(df) - nrow(var_es_dynamic$var01) + 1):nrow(df), ]
dates <- df$Date
df <- df[, -1]
factors <- colnames(df)

gather_series <- function(x, name = 'value') {
  result <- gather_(x, 'factor', name, factors)
  result$h <- rep(dates, 6)
  result$order <- factor(result$factor, factors)
  result
}

series <- gather_series(df)

ggplot(series, aes(x = h, y = value, group = factor)) +
  geom_line(aes(color = 'Realized Return')) +
  scale_y_continuous(labels = percent) +
  facet_grid(. ~ order)

ggplot(series, aes(x = h, group = factor)) +
  geom_line(aes(y = ES, color = '5% Dynamic ES'),
            data = gather_series(var_es_dynamic$var05, name = 'ES')) +
  geom_line(aes(y = ES, color = '5% Constant ES'),
            data = gather_series(var_es_constant$var05[1:287, ], name = 'ES')) +
  scale_y_continuous(labels = percent) +
  facet_grid(. ~ order)