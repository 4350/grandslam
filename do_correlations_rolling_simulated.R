#' Compute rolling 52 week correlations from copula
#' 


# Setup and pairs list ------------------------------------------------

rm(list = ls())

load('data/derived/weekly-estim.RData')

library(foreach)

MODEL_NAME <- 'std'
PATH <- file.path('data/derived/stdresid/dynamic', MODEL_NAME)

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

WINDOW = 52

# Functions ---------------------------------------------------------------

distribution <- lapply(list.files(PATH), function(p) {
  load(file.path(PATH, p))
  df <- data.frame(stdresid_t)
  colnames(df) <- c('Mkt.RF','HML','SMB','Mom','RMW','CMA')
  df
})

rolling_correlations_simulated <- function(distribution, window) {
  T_final <- length(distribution)
  
  foreach(t = window:T_final) %do% {
    # Subset the slices and stack them on top of each other
    subset <- distribution[(t - window + 1):t]
    subset <- bind_rows(subset)
    
    cor(subset)
  }
}


# Do ----------------------------------------------------------------------

simulated_roll_corr <- rolling_correlations_simulated(distribution, WINDOW)

# Save for plotting purposes ----------------------------------------------

roll_corr_simulated <- bind_rows(lapply(PAIRS, function(pair) {
  coef <- sapply(simulated_roll_corr, function(c) c[pair[1], pair[2]])
  data.frame(
    Date = tail(df.estim$Date, length(simulated_roll_corr)),
    factor1 = pair[1],
    factor2 = pair[2],
    model = MODEL_NAME,
    coef = coef
  )
}))

save(roll_corr_simulated, file = 'data/derived/rolling/roll_corr_simulated.RData')