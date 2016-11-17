#' Compute rolling 52 week correlations from copula
#' 
rm(list = ls())

library(foreach)

MODEL_NAME <- 'std'
PATH <- file.path('data/derived/stdresid/dynamic', MODEL_NAME)

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

simulated_roll_corr <- rolling_correlations_simulated(distribution, 52)
save(
  simulated_roll_corr,
  file = sprintf('data/correlations_rolling_simulated_%s.RData', MODEL_NAME)
)