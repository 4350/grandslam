#' Simulate from constant copula
#' 
#' Quickly written to optimize for this case. We generate stdresid
#' 

rm(list = ls())

simulate_constant <- function(model_name, n = 1e5) {
  # Load models
  load('data/derived/model_GARCH.RData')
  load(sprintf('data/derived/model_copula_%s.RData', model_name))
  
  # Name the copula model something simple
  var_name = sprintf('model.copula.%s', gsub('_', '.', model_name))
  model.Copula <- get(var_name)
  rm(list = var_name, var_name)
  
  # Get GARCH specs
  model.GARCH <- garch.fit2spec(model.GARCH)
  
  # Create output directory
  output_directory <- file.path('data/derived/stdresid', model_name)
  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
  
  # Simulation parameters
  Correlation <- model.Copula$Correlation[,, 1]
  uv <- dc.uv.dists(ncol(Correlation), model.Copula$params$dist.params)
  
  # Simulate some stuff
  shocks <- sim.c.rghyp(model.Copula, Correlation, n)
  stdresid <- shocks2stdresid(shocks, uv, model.GARCH)
  
  write.csv(
    stdresid,
    file.path(output_directory, sprintf('%d_1.csv', n))
  )
}

simulate_constant('constant_gauss')
simulate_constant('constant_ght')
simulate_constant('constant_ghskt')

