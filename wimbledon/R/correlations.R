#' Functionality related to rolling and threshold correlations
#' 

#' Threshold correlations
#' 
#' Takes a data frame, and returns estimates of threshold corr, lb, ub
#' using the factor names given
#' 
#' @param df Data frame of return series
#' @param factor1 String - factor number one
#' @param factor2 String - factor number two
#' 
#' @return result Matrix 91x3 of coef, lb, ub at percentiles 10-90
#' @export
correlations.threshold <- function(df, factor1, factor2) {
  
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

#' Rolling correlations
#' 
#' Takes a data frame and returns rolling 250-d correlations of pairs with lb ub,
#' using the factor names given
#' 
#' @param df Data frame of return series
#' @param factor1 String - factor number one
#' @param factor2 String - factor number two
#' 
#' @return result Matrix of coef, lb, ub for time series length less window size
#' @export
correlations.rolling <- function(df, factor1, factor2) {
  
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

#' Factor correlations apply
#' 
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
  correlations.factor.apply <- function(other.factors, value, statistic, df, thresholdtoggle) {
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
        return(rbind(correlations.threshold(df, value, other.factors)[statistic, ]))
      }
      return(rbind(correlations.rolling(df, value, other.factors)[statistic, ]))  
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

#' Plot rolling correlation. To be used in list apply
#' 
#' @param value Data frame for one factor of rolling correlation data
#'
#' @return gridded plot of rolling correlation graphs
#' @export
correlations.plot.rolling <- function(value) {
  ggplot(
    value,
    aes(x = Date, y = value, group = factor)
  ) + 
    geom_ribbon(aes(ymin = lb, ymax = ub),
                fill = "grey70") +
    geom_line(aes(color = factor)) +
    
    theme(legend.position="none") +
    ylab('') +
    coord_cartesian(ylim = c(-1, 1)) +
    #scale_x_date(date_breaks = "5 years", date_labels = "%y")+
    facet_grid(. ~ order)  
}

#' Plot threshold correlation. To be used in list apply
#' 
#' @param value Data frame for one factor of threshold correlation data
#'
#' @return gridded plot of threshold correlation graphs
#' @export
correlations.plot.threshold <- function(value) {
  ggplot(
    value,
    aes(x = qs, y = value, group = factor)
  ) + 
    geom_ribbon(aes(ymin = lb, ymax = ub),
                fill = "grey80") +
    geom_line(aes(color = factor)) +
    theme(legend.position="none") +
    ylab('') +
    coord_cartesian(xlim = c(0, 1), ylim = c(-0.5, 1)) +
    scale_x_continuous(labels = scales::percent) +
    facet_grid(. ~ order)  
}