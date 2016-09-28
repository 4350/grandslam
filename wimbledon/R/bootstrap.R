#' Generate stationary block Bootstrap sample indices
#'
#' @param B number of bootstrap samples 
#' @param T time series length per sample
#' @param length.block average block length
#'
#' @return TxB matrix of block of indices 1:T 
#' @export
bs.stationary.index <- function(B, T, length.block) {
  bs.block <- function(b) {
    # For each t, we start with a new block with 1/length.block probability
    new.block <- runif(T) < 1 / length.block
    new.block[1] <- TRUE
    
    i <- NA
    indices <- rep(NA, T)
    for (t in seq(T)) {
      if (new.block[t]) {
        # Pick random index to start block at
        i <- sample(seq(T), 1)
      }
      else {
        # Continue existing block; or loop back
        i <- ifelse(i == T, 1, i + 1)
      }
      
      indices[t] <- i
    }
    
    indices
  }
  
  sapply(seq(B), bs.block)
}