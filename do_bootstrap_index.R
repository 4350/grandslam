#' Bootstrap Block Selection
#' 
#' The only thing this code does is to select blocks for the bootstrap. These
#' are later used by the actual bootstrapping procedure.

# Setup ------------------------------------------------------------------

rm(list = ls())

set.seed(403)

BLOCK_LENGTH <- 104
BS_INDEX <- bs.stationary.index(1000, 2766, BLOCK_LENGTH)

save(BS_INDEX, file = 'data/derived/bootstrap/bs_index.RData')
