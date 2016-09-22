".besselM3" <- function(lambda = 9/2, x = 2, logvalue = FALSE) {
  if(all(abs(lambda) == 0.5)){
    ## Simplified expression in case lambda == 0.5
    if (!logvalue){
      res <- sqrt(pi/(2 * x)) * exp(-x)
    }else{
      res <- 0.5 * log(pi/(2 * x)) - x
    }
  }else{
    if (!logvalue){
      res <- besselK(x, lambda)
    }else{
      res <- log( besselK(x, lambda, expon.scaled = TRUE) ) - x
    }
  }
  return(res)
}

#' Special multivariate density function for asymmetric t-distribution
#' 
#' Calling this function directly instead of going through dghyp shaves about
#' 500 milliseconds per run (for big T). It skips most error checking, assumes
#' mu = 0
#' 
#' Code copied from internal dghypmv and stripped of essentially all error
#' checking.
dghsst <- function(x, nu, gamma, sigma) {
  lambda = -nu / 2
  chi = nu - 2
  
  d <- dim(x)[2]
  n <- dim(x)[1]
  
  det.sigma <- det(sigma)
  inv.sigma <- solve(sigma)
  #Q <- mahalanobis(x, FALSE, inv.sigma, inverted = TRUE)
  Q <- rowSums(x %*% inv.sigma * x)
  
  g <- inv.sigma %*% gamma 
  skewness.scaled <- as.vector(x %*% g)
  skewness.norm <- t(gamma) %*% g
  
  lambda.min.d.2 <- lambda - d / 2
  interm <- sqrt((chi + Q) * skewness.norm)
  
  log.const.top <- -lambda * log(chi) - lambda.min.d.2 * log(skewness.norm)
  log.const.bottom <- d / 2 * log(2 * pi) + 0.5 * log(det.sigma) +
    lgamma(-lambda) - (lambda + 1) * log(2)
  log.top <-
    .besselM3(lambda.min.d.2, interm, logvalue = TRUE) +
    skewness.scaled
  log.bottom <- -lambda.min.d.2 * log(interm)
  
  log.const.top + log.top - log.const.bottom - log.bottom
}