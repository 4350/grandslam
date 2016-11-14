#' Give a function that computes CDB for period t
#' 
#' Call the returned function with portfolio weights to get the CDB
#' for this period. Perfect for optimization.
#'
#' @param q The significance level
#' @param returns Returns for assets in period t
#'
#' @return Negative CDB
#' @export
cdb_fn <- function(q, returns) {
  # Compute the risk measures of each factor separately
  var <- apply(returns, 2, function(r) -quantile(r, q))
  es <- lapply(seq_along(var), function(i) {
    r <- returns[, i]
    -mean(r[r <= -var[i]])
  })
  es <- unlist(es)
  
  function(weights) {
    portfolio_returns <- returns %*% weights
    portfolio_var <- -quantile(portfolio_returns, q, names = F)
    portfolio_es <- -mean(portfolio_returns[portfolio_returns <= -portfolio_var])
    
    es_ub <- sum(weights * es)
    es_lb <- portfolio_var
    ((es_ub - portfolio_es) / (es_ub - es_lb))
  }
}