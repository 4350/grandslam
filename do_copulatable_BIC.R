
# Startup -----------------------------------------------------------------
rm(list = ls())
library(devtools)
load_all('wimbledon')
load_all('australian')

# Constant copulas --------------------------------------------------------

# List of data frames
load('data/derived/copula/full_constant.RData')
load('data/derived/copula/full_dynamic.RData')

copulatable <- function(copula_list) {
  out.data <- lapply(copula_list, function(model) {
    N = ncol(model$fit@dynamics@Omega)
    
    # NOTE HARDCODED
    ################
    T = 2766
    
    Parameters = model$optimized$par
    Observations = T
    ll = model$ll
    nParams = length(Parameters) + (N * (N-1) / 2)
    BIC = -2 * ll + nParams * log(T)
    Persistence = tryCatch(model$optimized$par['alphabeta.alpha'], error = function(err) NA) + 
      tryCatch(model$optimized$par['alphabeta.beta'], error = function(err) NA)
    
    out.data <- data.frame(value =
                             unlist(
                               list(
                                 Parameters,
                                 Observations = Observations,
                                 ll = ll,
                                 nParams = nParams,
                                 BIC = BIC,
                                 Persistence = Persistence
                               )
                             )
    )
    
    out.data <- round(out.data, digits = 3)
  })
  
  print(out.data, digits = 3)  
}

copulatable(constant_copula_fit)
copulatable(dynamic_copula_fit)
