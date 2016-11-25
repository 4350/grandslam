#' Takes set of weights and returns some portfolio results
#' using realized data
#' 
#' @param weights matrix of weights, should be named
#' @param distribution full dist of returns
#' @param q CDB percentage VaR cutoff
#' @param selectors vector of strings with facotrs included
#'
#' @return list of portfolio results
#' @export
portfolio_metrics <- function(weights, distribution, q, selectors) {
  
  load('data/derived/weekly-estim.RData')
  
  T = dim(distribution)[3]

  # Use realized returns and dates
  dates <- tail(df.estim[,'Date'], T)
  # Subset realized returns using selectors
  realized <- tail(df.estim[, selectors], T)
  
  # Calculate the portfolios realized return
  portfolio_return <- rowSums(weights * realized)
  
  # Calculate the CDB and VaR
  cdb_var_es_results <- cdb_var_es(q, distribution, weights)
  
  # List of results
  results <- c(
    list(Date = dates$Date,
         portfolio_return = portfolio_return),
    cdb_var_es_results
  )
  
}