#' Functionality related to rolling correlations
#'

#' Rolling correlations
#'
#' Takes a data frame, and returns estimates of rolling corr, lb, ub
#' using the factor names given
#'
#' 
#' @param pair Two strings in character vector
#' @param window length of roll window
#' @param df Data frame of return series
#' 
#' @return result GGplot-tidy data frame w coef, lb, ub, unc_coef, factors
#' @export
correlations_rolling <- function(pair, window, df) {
  
  factor1 <- pair[[1]]
  factor2 <- pair[[2]]
  
  # Cut df so that it is as long as simulations from dist
  df <- df[2:nrow(df),]
  dates <- df[,'Date']
  df <- df[,pair]
  
  T_final <- nrow(df)
  
  result <- foreach(t = window:T_final) %do% {
    # Subset the data and get correlation
    subset <- df[(t - window + 1):t, pair]
    
    result <- cor.test(subset[,factor1], subset[,factor2])
    
    # Return coef, lb, ub, unc_coef
    result <- data.frame(coef = result$estimate[[1]],
                         lb = result$conf.int[1],
                         ub = result$conf.int[2],
                         unc_coef = cor(df[,factor1], df[,factor2])[[1]]
                         )
    
  }
  
  
  result <- data.frame(
    bind_rows(result),
    Date = dates[window:T_final],
    factor1 = factor1,
    factor2 = factor2
  )
  
  return(result)
}

#' Get rolling correlation data frame for plotting for a
#' set of factor pairs (in pairs list)
#' @param pairs list of vectors with the two factor strings
#' @param window length of roll window
#' @param df to run the threshold correlation on
#' @return results rows bound together for all pairs, tidy data
#' @export
rolling_correlations_pairs <- function(pairs, window, df) {
  # List apply pairs in threshold correlation function
  results <- lapply(pairs, function(pair, window, df) {
    correlations_rolling(pair, window, df)
  },
  window = window, 
  df = df
  )
  
  bind_rows(results)
  
}