#' Function to get values for Standardized Residuals Empirical Denisty graph
#'
#' @param x fitted model object
#' @return out.df Data frame with histogram bins, limits, normal and sstd lines
#' 
#' @export
#' 
empirical.density = function(x, ...)
{
  vmodel  = x@model$modeldesc$vmodel
  zseries = as.numeric(residuals(x, standardize=TRUE))
  distribution = x@model$modeldesc$distribution
  idx = x@model$pidx
  pars  = x@fit$ipars[,1]
  skew  = pars[idx["skew",1]]
  shape = pars[idx["shape",1]]
  if(distribution == "ghst") ghlambda = -shape/2 else ghlambda = pars[idx["ghlambda",1]]	
  xlim 	= c(min(zseries), max(zseries))
  result 	= hist(x = zseries, col = "grey", border = "white",
                 breaks = "Scott", main = "Empirical Density of Standardized Residuals", xlim = xlim, ylim = c(0,0.6),
                 probability = TRUE, ylab="Probability", cex.main = 0.8, 
                 plot = FALSE, ...)
  s = result$mids
  y	= ddist(distribution, s, lambda = ghlambda, skew = skew, shape = shape)
  normline = dnorm(s,0,1)
  histdensity = result$density
  
  # Output data frame with midpoints of bins (x), normal and sstd line and densities
  out.df <- gather(
    data.frame(x = s, Normal = normline, Sstd = y, empirical = histdensity),
    'Distribution',
    'values',
    2:4
  )
  
  
  
}