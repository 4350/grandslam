#' Takes a data frame and returns rolling 250-d correlations of pairs with lb ub,
#' using the factor names given
#' 
#' @param df Data frame of return series
#' @param factor1 String - factor number one
#' @param factor2 String - factor number two
#' 
#' @return result Matrix of coef, lb, ub for time series length less window size
#' @export
rollingCorr <- function(df, factor1, factor2) {
  
  corrvalues <- df %>%
    .[,c(factor1, factor2)] %>%
    as.matrix() %>%
    
    rollapply(250, function(x) {
      result <- cor.test(x[, 1], x[, 2], method = 'pearson')
      corresult <- c(result$estimate,
                     result$conf.int[1],
                     result$conf.int[2])
    }, 
    by.column=FALSE, align="right")
  
  corrvalues <- t(corrvalues)
  rownames(corrvalues) <- c('coef', 'lb', 'ub')
  return(corrvalues)
  
}