#' Takes a data frame, and returns estimates of threshold corr, lb, ub
#' using the factor names given
#' 
#' @param df Data frame of return series
#' @param factor1 String - factor number one
#' @param factor2 String - factor number two
#' 
#' @return result Matrix 91x3 of coef, lb, ub at percentiles 10-90
#' @export
thresholdCorrelations <- function(df, factor1, factor2) {
  
  qs <- seq(0.10, 0.90, by=0.01)
  
  result <- sapply(qs, function(q, df, factor1, factor2) {
    # Extract vectors
    x = simplify2array(df[, factor1])
    y = simplify2array(df[, factor2])
    # Find empirical percentiles
    quantile.x <- quantile(x, q)
    quantile.y <- quantile(y, q)

    # Check which of x and y to include in correlation estimation
    if (q < 0.50) {
      inc <- x < quantile.x & y < quantile.y
    }
    else {
      inc <- x >= quantile.x & y >= quantile.y
    }

    # Estimate correlation
    result <- cor.test(x[inc], y[inc])

    # Return coef, lb, ub
    corresult <- c(result$estimate,
                   result$conf.int[1],
                   result$conf.int[2])

    return(corresult)
  },
  df = df, 
  factor1 = factor1,
  factor2 = factor2
  )
  
  rownames(result) <- c('coef', 'lb', 'ub')
  return(result)
}