
# Setup ------------------------------------------------------------------

rm(list = ls())

library(dplyr)
library(tictoc)
library(devtools)
load_all('wimbledon')

load('data/derived/garch/model_GARCH_chosen_stdres.RData')
load('data/derived/weekly-estim.RData')
df.estim <- as.data.frame(df.estim)

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

WINDOW <- 52

# Do and save stdres data set ---------------------------------------------

roll_corr_stdres <- rolling_correlations_pairs(PAIRS, WINDOW, df.stdres)
roll_corr_stdres$model <- 'stdres'

save(roll_corr_stdres, file = 'data/derived/rolling/roll_corr_stdres.RData')


# Do and save returns data set ---------------------------------------------

roll_corr_returns <- rolling_correlations_pairs(PAIRS, WINDOW, df.estim)
roll_corr_returns$model <- 'returns'

save(roll_corr_returns, file = 'data/derived/rolling/roll_corr_returns.RData')


