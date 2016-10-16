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
#' @return result Matrix 81x3 of coef, lb, ub at percentiles 10-90
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
#' Takes a data frame and returns rolling n-period correlations of pairs with lb ub,
#' using the factor names given
#'
#' @param df Data frame of return series
#' @param factor1 String - factor number one
#' @param factor2 String - factor number two
#' @param window Length of correlation roll
#' 
#' @return result Matrix of coef, lb, ub for time series length less window size
#' @export
correlations.rolling <- function(df, factor1, factor2, window) {
  
  corrvalues <- df %>%
    .[,c(factor1, factor2)] %>%
    as.matrix() %>%
    
    rollapply(width = window, function(x) {
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
#' @param window Length of correlation window in rolling
#' 
#' @return result Data frame (tidy) with factor and value of statistic
#' @export
  correlations.factor.apply <- function(other.factors, value, statistic, df, thresholdtoggle, window = NULL) {
  result = sapply(
    other.factors,
    function(other.factors, value, statistic, df, thresholdtoggle, window) {
      # Check if same column for value factor as for other factor, then return NA
      if (other.factors == value) {
        # Check if in rolling or threshold
        if (thresholdtoggle == 1) {
          return(matrix(NA, 1, 81))
        }
        return(matrix(NA, 1, dim(df)[1]-window+1))
      }
      # Else return the result depending on which statistic is requested
      # Check if in rolling or threshold
      if (thresholdtoggle == 1) {
        return(rbind(correlations.threshold(df, value, other.factors)[statistic, ]))
      }
      return(rbind(correlations.rolling(df, value, other.factors, window)[statistic, ]))  
    },
    value = value,
    statistic = statistic,
    df = df,
    thresholdtoggle = thresholdtoggle,
    window = window
  )

  # Gather result in tidy fashion
  result <- gather(
    data.frame(result),
    'factor',
    'value',
    1:length(other.factors)
  )
  # Add facetvalue, one of the value factors, to do even better ggplot
  factors.value = c("HML", "RMW", "CMA")
  factors.all   = c("Mkt.RF", "HML", "SMB", "Mom", "RMW", "CMA")
  
  result$facetvalue <- value
  result$order2 <- factor(result$facetvalue, factors.value)
  return(result)
}

#' Calculate list of threshold correlations
#'
#' @param df Data frame (TxBootruns)xN of simulated data or empirical data (shorter)
#' @param simulatetoggle Dummy 1 if simulated series, else considered empirical
#'
#' @return out.list a list of three (one for each value factor)
#'  data frames with threshold correlations for each factor,
#'  also holds factor ordering for ggplot
#' @export
th_corr <- function(df, simulatetoggle) {
  # Get some size of data
  N <- ncol(df)
  T <- nrow(df)

  # Set up factor groups ----
  factors.value <- c("HML", "RMW", "CMA")
  factors.all <- c("Mkt.RF", "HML", "SMB", "Mom", "RMW", "CMA")

  # Create list to hold simulated corrs
  out.list = list()

  # Populate output matrix with point estimates for all value factors vs all factors
  for (value in factors.value) {
    # Get point estimates for each pair value <-> other factor
    # If simulated data set
    if (simulatetoggle == 1) {
      out.df <- data.frame(
        simcoef = correlations.factor.apply(factors.all, value, 'coef', df, 1)[,'value']
      )
    } else {
    # Else empirical data set
      out.df <- data.frame(
        qs = seq(0.10,0.90,0.01),
        correlations.factor.apply(factors.all, value, 'coef', df, 1),
        lb = correlations.factor.apply(factors.all, value, 'lb', df, 1)[,'value'],
        ub = correlations.factor.apply(factors.all, value, 'ub', df, 1)[,'value']
      )
      # Ordering for ggplot
      out.df$order <- factor(out.df$factor, levels = factors.all)
    }
    # Put in list
    out.list[[value]] <- out.df
  }
  return(out.list)
}

#' Calculate list of rolling n-period correlations
#' 
#' @param df Data frame of empirical weekly data TxN
#' @param df.date Optional, for date numbering
#' @param window Length of corr window
#'
#' @return out.list a list of three (one for each value factor)
#'  data frames with rolling correlations for each factor,
#'  also holds factor ordering for ggplot
#' @export
roll_corr <- function(df, df.date = NULL, window) {
  # Get some size of data
  N <- ncol(df)
  T <- nrow(df)
  nCorrs <- dim(df.estim)[1]-window+1
  
  # Check if dates are specified or whether to use numbers
  if(is.null(df.date)) {
    df.date = seq(nCorrs)
  }
  
  # Set up factor groups ----
  factors.value <- c("HML", "RMW", "CMA")
  factors.all <- c("Mkt.RF", "HML", "SMB", "Mom", "RMW", "CMA")
  
  # Create list to hold corrs
  out.list = list()
  
  # Populate output matrix with point estimates for all value factors vs all factors
  for (value in factors.value) {
    # Get point estimates for each pair value <-> other factor
    out.df <- data.frame(
      Date = tail(df.date, nCorrs),
      correlations.factor.apply(factors.all, value, 'coef', df, 0, window),
      lb = correlations.factor.apply(factors.all, value, 'lb', df, 0, window)[,'value'],
      ub = correlations.factor.apply(factors.all, value, 'ub', df, 0, window)[,'value']
    )
    
    # Ordering for ggplot
    out.df$order <- factor(out.df$factor, levels = factors.all)
    
    # Put in list
    out.list[[value]] <- out.df
  }
  return(out.list)
}

#' Do threshold correlation from simulated data directory
#' Consolidates all simulation runs and standardizes the returns
#' 
#' @param file The directory of simulated data for this distribution
#' 
#' @return list of dataframes with threshold correlations
#' @export
do.th.sim <- function(file) {
  # Set up factor groups
  factors.value = c("HML", "RMW", "CMA")
  factors.all   = c("Mkt.RF", "HML", "SMB", "Mom", "RMW", "CMA")
  
  df <- get.simulation.results(file)
  df <- df[, -1]
  colnames(df) <- factors.all
  
  th_corr(df, 1)
}

#' Consolidates all simulation runs into one matrix
#' 
#' @param dir the directory where simulation csvs are
#' 
#' @return df consolidated df of all simulations runs' returns, sigma, residuals
#' @export
get.simulation.results <- function(dir) {
  filenames <- file.path(dir, list.files(dir))
  do.call('rbind', lapply(filenames, read.csv))
}