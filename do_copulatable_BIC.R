
# Startup -----------------------------------------------------------------
rm(list = ls())
library(devtools)
load_all('wimbledon')


# Constant copulas --------------------------------------------------------

# List of data frames
load('data/derived/model_copula_constant_gauss.RData')
load('data/derived/model_copula_constant_ght.RData')
load('data/derived/model_copula_constant_ghskt.RData')
load('data/derived/model_copula_dynamic_gauss.RData')
load('data/derived/model_copula_dynamic_ght.RData')
load('data/derived/model_copula_dynamic_ghskt.RData')


model_list <- list(constant_gauss = model.copula.constant.gauss, 
                   constant_ght = model.copula.constant.ght, 
                   constant_ghskt = model.copula.constant.ghskt,
                   dynamic_gauss = model.copula.dynamic.gauss,
                   dynamic_ght = model.copula.dynamic.ght,
                   dynamic_ghskt = model.copula.dynamic.ghskt
)

out.data <- lapply(model_list, function(model) {
  N = ncol(model$shocks)
  T = nrow(model$shocks)
  
  Parameters = model$params
  Observations = T
  ll = model$ll
  nParams = length(unlist(model$params)) - (model$params$alpha == 0) - (model$params$beta == 0) + (N * (N-1) / 2)
  BIC = -2 * ll + nParams * log(T)
  Persistence = model$params$alpha + model$params$beta
  
    
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