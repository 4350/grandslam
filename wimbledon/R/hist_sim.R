#'histsim functionality
#'

#' Function to estimate GARCH models with full data set, leaving some data OOS
#' @param df data frame of returns to estimate garch columnwise on
#' @param H number of days to be left out-of-sample
#' 
#' @return garch list of garch fits objects, one for each column
#' @export
estimate.garch <- function(df, H) {
  garch <- lapply(names(kGARCHModels), function(name, H) {
    model <- kGARCHModels[[name]]
    data <- as.data.frame(df[, name])
    ugarchfit(model, data, out.sample = H, solver = 'hybrid')
  },
  H = H)
  names(garch) <- names(kGARCHModels)
  
  garch
}

#' fc_List makes a list of forecast data frames
#' @param garch.fits list of garch fit objects with OOS data of length H
#' @param H horizon of forecast
#' 
#' @return out.list a list containing forecast df for each ret series
#' @export
fc_List <- function(garch.fits, H) {
  out.list = list()
  n.roll = H - 1
  out.list <- lapply(garch.fits,
                     function(fFit, n.roll) {
                       ugarchforecast(
                         fitORspec = fFit,
                         n.ahead = 1,
                         n.roll = n.roll
                       )
                     },
                     n.roll = n.roll
  )
}

#' Function to calculate volatility updated returns
#' see Hull & White (1998)
#' 
#' @param h Index in the forecast period, h<=H
#' @param q Lower tail probability of interest
#' @param df Data frame of returns (full sample)
#' @param sigma Consolidated matrix of conditional volatilities (T+H)*N, both in and out of sample
#' @param T End of estimation period
#' 
#' @return var.star Vector length N. 1-step ahead VaR for forecast index h
#' @export
volatilityupdate <- function(h, q, df, sigma, T) {
  
  # Standard returns up to last period
  std.ret <- df[1:(T+h-1), ] / sigma[1:(T+h-1), ]
  r.star <- apply(std.ret, 1, function(r, h) {
    # Multiply each time t standardized rets by next forecast of vol
    r * sigma[(T+h),]
  },
  h = h)
  
  # Turn it right
  r.star <- t(r.star)
}

#' Function to calculate VaR using Hull & White vol updated returns
#' 
#' @param h Index in the forecast period, h<=H
#' @param q Lower tail probability of interest
#' @param df Data frame of returns (full sample)
#' @param sigma Consolidated matrix of conditional volatilities (T+H)*N, both in and out of sample
#' @param T End of estimation period
#' 
#' @return var.star
#' @export
var.HW <- function(h, q, df, sigma, T) {
  r.star <- volatilityupdate(h, q, df, sigma, T)
  #Calculate VaRs
  var.star <- apply(r.star, 2, function(r, q) {
    quantile(r, probs = q)
  },
  q = q
  )
  
  var.star
}

#' Function to calculate ES using Hull & White vol updated returns
#' 
#' @param h Index in the forecast period, h<=H
#' @param q Lower tail probability of interest
#' @param df Data frame of returns (full sample)
#' @param sigma Consolidated matrix of conditional volatilities (T+H)*N, both in and out of sample
#' @param T End of estimation period
#' 
#' @return es.star
#' @export
es.HW <- function(h, q, df, sigma, T) {
  r.star <- volatilityupdate(h, q, df, sigma, T)
  #Calculate VaRs
  var.star <- apply(r.star, 2, function(r, q) {
    quantile(r, probs = q)
  },
  q = q
  )
  # Count rows and series
  kRows <- nrow(r.star)
  N <- ncol(df)
  # This matrix can be compared to r.star 
  var.star.matrix <- matrix(rep(var.star, kRows), kRows, N, byrow = TRUE)
  # Logical indexing over columns
  es.star <- sapply(seq(N), function(c) {
    mean(r.star[(r.star[,c]<var.star.matrix[,c]),c])
  })
  
}



#' Historical simulation of VaR using volatility updated returns
#' 
#' @param q Lower CDF probability of interest
#' @param std.ret Consolidated matrix of standardized returns (T+H)xN, both in and out of sample
#' @param sigma Consolidated matrix of conditional volatilities (T+H)*N, both in and out of sample
#' @param T End of estimation period
#' 
#' @return hs.var df of historically simulated VaR
#' @export
hist.sim.var <- function(q, df, sigma, T) {
  H = nrow(df) - T
  hs.var = sapply(seq(1:H), function(h, q, df, sigma, T) {
    var.HW(h, q, df, sigma, T)
  },
  q = q,
  df = df,
  sigma = sigma, 
  T = T
  )
  hs.var = data.frame(t(hs.var))
}

#' Historical simulation of ES using volatility updated returns
#' 
#' @param q Lower CDF probability of interest
#' @param std.ret Consolidated matrix of standardized returns (T+H)xN, both in and out of sample
#' @param sigma Consolidated matrix of conditional volatilities (T+H)*N, both in and out of sample
#' @param T End of estimation period
#' 
#' @return hs.es df of historically simulated ES
#' @export
hist.sim.es <- function(q, df, sigma, T) {
  H = nrow(df) - T
  hs.es = sapply(seq(1:H), function(h, q, df, sigma, T) {
    es.HW(h, q, df, sigma, T)
  },
  q = q,
  df = df,
  sigma = sigma, 
  T = T
  )
  hs.es = data.frame(t(hs.es))
}

#' Get tidy VaR df for plotting purposes
#' @param q Lower tail probability of interest
#' @param df Data frame of returns (full sample)
#' @param sigma In-sample and out-of-sample conditional vol forecasts
#' @param T time when the in-sample ends and out-of-sample begins
#' @param kGARCHModels list of fixed models used, for naming purposes
#' 
#' @return out.df tidy VaR data frame
#' @export
hist.sim.tidy.var <- function(q, df, sigma, T, kGARCHModels) {
  H = nrow(df) - T
  N = ncol(df)
  # Run simulation
  out.df <- hist.sim.var(q, df, sigma, T)
  # Long format data
  out.df <- gather(out.df, 'factor', 'VaR', 1:N)
  # Add horizon "x values"
  out.df$h <- rep(seq(1,H), N)
  # Order data set for plot
  out.df$order <- factor(out.df$factor, names(kGARCHModels))
  out.df
}

#' Get tidy ES df for plotting purposes
#' @param q Lower tail probability of interest
#' @param df Data frame of returns (full sample)
#' @param sigma In-sample and out-of-sample conditional vol forecasts
#' @param T time when the in-sample ends and out-of-sample begins
#' @param kGARCHModels list of fixed models used, for naming purposes
#' 
#' @return out.df tidy ES data frame
#' @export
hist.sim.tidy.es <- function(q, df, sigma, T, kGARCHModels) {
  H = nrow(df) - T
  N = ncol(df)
  # Run simulation
  out.df <- hist.sim.es(q, df, sigma, T)
  # Name the columns as models
  colnames(out.df) <- names(kGARCHModels)
  # Long format data
  out.df <- gather(out.df, 'factor', 'ES', 1:N)
  # Add horizon "x values"
  out.df$h <- rep(seq(1,H), N)
  # Order data set for plot
  out.df$order <- factor(out.df$factor, names(kGARCHModels))
  out.df
}

#' Kupiec test and Christoffersen test for HSHW VaR. Saves a csv table.
#' 
#' @param q Lower tail probaiblity of interest
#' @param df Data frame (full sample)
#' @param sigma In-sample and out-of-sample conditional vol
#' @param T End of in-sample period
#' @param kGARCHModels List of fixed garch models (for naming)
#' 
#' @return out.matrix 12xN matrix of test statistics (for each return series)
#' @export
kc_test.var <- function(q, df, sigma, T, kGARCHModels) {
  N = ncol(df)
  df.var <- hist.sim.var(q, df, sigma, T)
  df.emp <- as.data.frame(df[(T+1):nrow(df), ])
  
  out.matrix <- matrix(NA, 12, N)
  colnames(out.matrix) <- names(kGARCHModels)
  
  for (c in seq(N)) {
    result <- VaRTest(alpha = q, df.emp[,c], df.var[,c])
    out.matrix[,c] <- simplify2array(result[1:12])
  }
  
  rownames(out.matrix) <- names(result[1:12])  
  write.table(out.matrix, file = paste('output/VaR/kupiec', 'T', T, 'q', q, 'csv', sep = '.'), sep = ',')
  return(out.matrix)
}

#' Get ES
#' @param q Lower tail probability of interest
#' @param df Data frame of returns (full sample)
#' @param sigma In-sample and out-of-sample conditional vol forecasts
#' @param T time when the in-sample ends and out-of-sample begins
#' 
#' @return 
#' @export
function(q, df, sigma, T) {
  H = nrow(df) - T
  N = ncol(df)
  # Get VaR matrix
  df.var <- hist.sim.var(q, df, sigma, T)
  # Get return data for same period
  df.emp <- df[(T+1):nrow(df),]
  # Compare return data to VaR to see where return is in the tail
  sapply(1:6, function(c) {
    mean(df.emp[(df.emp[,c]<df.var[,c]),c])
  })
  
  
  
}

