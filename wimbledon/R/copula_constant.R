#' Build a correlation matrix from a list of correlations
#' 
#' The order of correlations shouldn't matter, but they become are listed
#' column by column in the lower triangle
#'
#' @param correlations 
#'
#' @return NxN correlation matrix
#' @export
as.correlation.matrix <- function(correlations) {
  # Solve N * (N - 1) / 2 = length(correlations) for N
  N <- 0.5 * (sqrt(4 * 2 * length(correlations) + 1) + 1)
  
  c <- diag(1, N)
  c[lower.tri(c)] <- correlations
  c[upper.tri(c)] <- t(c)[upper.tri(c)]
  c
}

cc.ll.total <- function(u, mvdist) {
  # Extract the univariate distributions from the multivariate object
  uvdists <- lapply(seq(ncol(u)), function(i) {
    # sigma is always unit diagonal, but let's do it the right way
    sigma <- sqrt(diag(mvdist@sigma)[i])
    
    if (ghyp:::.is.gaussian(mvdist)) {
      return(ghyp::gauss(sigma = sigma))
    }
    else {
      return(ghyp::student.t(
        nu = mvdist@chi + 2,
        gamma = mvdist@gamma[i],
        sigma = sigma
      ))
    }
  })
  
  # Inverse the uniforms into their copula shocks
  x <- sapply(seq(ncol(u)), function(i) {
    ghyp::qghyp(u[, i], uvdists[[i]], method = 'splines')
  })
  
  ll.joint <- ghyp::dghyp(x, mvdist, logvalue = T)
  ll.marginal <- rowSums(sapply(seq(ncol(u)), function(i) {
    ghyp::dghyp(x[, i], uvdists[[i]], logvalue = T)
  }))
  
  # Sklar
  sum(ll.joint) - sum(ll.marginal)
}