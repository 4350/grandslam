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

cdb_var_es <- function(q, distribution, weights) {
  
  times = nrow(weights)
  cdb <- c(NA, T)
  var <- c(NA, T)
  es <- c(NA, T)
  
  for(t in 1:times) {
    fn <- cdb_fn(q, distribution[,,t])
    cdb[t] <- fn(weights[t,])
    fn_var <- var_fn(q, distribution[,,t])
    var[t] <- fn_var(weights[t,])
    fn_es <- es_fn(q, distribution[,,t])
    es[t] <- fn_es(weights[t,])
  }
  
  list(cdb = cdb, var = var, es = es)
}

var_fn <- function(q, returns) {
  # Compute the risk measures of each factor separately
  var <- apply(returns, 2, function(r) -quantile(r, q))
  es <- lapply(seq_along(var), function(i) {
    r <- returns[, i]
    -mean(r[r <= -var[i]])
  })
  es <- unlist(es)
  
  function(weights) {
    portfolio_returns <- returns %*% weights
    -quantile(portfolio_returns, q, names = F)
  }
}

es_fn <- function(q, returns) {
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
    -mean(portfolio_returns[portfolio_returns <= -portfolio_var])
  }
}