
# Setup ------------------------------------------------------------------

rm(list = ls())

library(dplyr)
library(tictoc)
library(devtools)
load_all('wimbledon')

load('data/derived/weekly-estim.RData')
df.estim <- as.data.frame(df.estim)
load('data/derived/garch/model_GARCH_chosen_stdres.RData')

# Setup function arguments ------------------------------------------------

PAIRS <- list(
  c('Mkt.RF', 'HML'),
  c('Mkt.RF', 'CMA'),
  c('Mkt.RF', 'RMW'),
  c('Mom', 'HML'),
  c('Mom', 'CMA'),
  c('Mom', 'RMW'),
  c('SMB', 'HML'),
  c('SMB', 'CMA'),
  c('SMB', 'RMW'),
  c('HML', 'CMA'),
  c('HML', 'RMW'),
  c('CMA', 'RMW')
)

QS <- seq(0.10, 0.90, 0.01)


# Do and save stdres data set ---------------------------------------------

th_corr_stdres <- threshold_correlations_pairs(PAIRS, QS, df.stdres)
th_corr_stdres$model <- 'stdres'

save(th_corr_stdres, file = 'data/derived/threshold/th_corr_stdres.RData')


# Do and save returns data set ---------------------------------------------

th_corr_returns <- threshold_correlations_pairs(PAIRS, QS, df.estim)
th_corr_returns$model <- 'returns'

save(th_corr_returns, file = 'data/derived/threshold/th_corr_returns.RData')


