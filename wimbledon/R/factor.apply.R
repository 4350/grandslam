#' Function to apply other factors on a value factor and return a statistic of the correlation
#' Used in for loop for each value factor and each statistic (coef, lb, ub)
#' 
#' @param other.factors List of all factors
#' @param value String, value factor used
#' @param statistic String, indicating the required statistic to get
#' @param df data frame used
#' @param threshold 1 for threshold, 0 for rolling
#' 
#' @return result Data frame (tidy) with factor and value of statistic
#' @export
factor.apply <- function(other.factors, value, statistic, df, thresholdtoggle) {
  result = sapply(
    other.factors,
    function(other.factors, value, statistic, df, thresholdtoggle) {
      # Check if same column for value factor as for other factor, then return NA
      if (other.factors == value) {
        # Check if in rolling or threshold
        if (thresholdtoggle == 1) {
          return(matrix(NA, 1, 81))
        }
          return(matrix(NA, 1, dim(df)[1]-250+1))
      }
      # Else return the result depending on which statistic is requested
      # Check if in rolling or threshold
      if (thresholdtoggle == 1) {
        return(rbind(thresholdCorrelations(df, value, other.factors)[statistic, ]))
      }
        return(rbind(rollingCorr(df, value, other.factors)[statistic, ]))  
    },
    value = value,
    statistic = statistic,
    df = df,
    thresholdtoggle = thresholdtoggle
  )

  # Gather result in tidy fashion
  result <- gather(
    data.frame(result),
    'factor',
    'value',
    1:length(factors.all)
  )
}