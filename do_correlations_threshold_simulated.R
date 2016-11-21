# Setup ------------------------------------------------------------------

rm(list = ls())

library(dplyr)
library(tictoc)
library(devtools)
load_all('wimbledon')



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

MODELS <- list(
  'norm',
  'std',
  'ghst'
)

QS <- seq(0.10, 0.90, 0.01)

DYNAMIC <- 'constant'


# Do and save data set ----------------------------------------------------

tic(sprintf(
  'Threshold correlations for %d pairs, %d models, between %.2f and %.2f quantiles',
  length(PAIRS),
  length(MODELS),
  QS[1],
  tail(QS,1)
  )
)
th_corr_simulated <- threshold_correlations_models(DYNAMIC, MODELS, PAIRS, QS)
toc()

save(th_corr_simulated, file = sprintf('data/derived/threshold/th_corr_%s_simulated.RData',
                            DYNAMIC)
     )
