#' Compute rolling 52 week correlations from copula
#' 
rm(list = ls())

library(foreach)

MODEL_NAME <- 'std'
load(sprintf('data/derived/distributions/full_dynamic_%s_10000.RData', MODEL_NAME))

rolling_correlations_simulated <- function(distribution, window) {
  T_final <- dim(distribution)[3]
  
  foreach(t = window:T_final) %do% {
    # Subset the slices and stack them on top of each other
    subset <- distribution[,, (t - window + 1):t]
    subset <- bind_rows(
      lapply(seq(dim(subset)[3]), function(i) data.frame(subset[,, i]))
    )
    
    cor(subset)
  }
}

simulated_roll_corr <- rolling_correlations_simulated(distribution, 52)
save(
  simulated_roll_corr,
  file = sprintf('data/correlations_rolling_simulated_%s.RData', MODEL_NAME)
)