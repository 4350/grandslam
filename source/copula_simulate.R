garch_path_t <- function(t, filtered_copula) {
  # Get simulated stdresid for this ("next") period
  stdresid_t <- .simulate_stdresid_t(t, filtered_copula)
  
  foreach (i = seq_along(model.GARCH), .combine = 'cbind') %do% {
    path <- garch.path.it(i, t, specs, filtered_series, stdresid_t)
    t(path@path$seriesSim)
  }
}

#' Simulate standardized residuals from copula
#' 
#' Uses a memoised version of the underlying simulator if using a constant
#' copula.
.simulate_stdresid_t <- function(t, filtered_copula) {
  if (copula_is_constant(filtered_copula$spec)) {
    return(.do_simulate_stdresid_t_constant(filtered_copula))
  }
  
  .do_simulate_stdresid_t(t, filtered_copula)
}

#' Simulate t+1 uniforms from copula
.copula_simulate_t <- function(t, filtered_copula) {
  Q_t <- filtered_copula$Q[,, t]
  shocks_t <- filtered_copula$shocks[t, ]
  
  t(simplify2array(
    copula_simulate(filtered_copula$spec, 1, N_RANDOM, Q_t, shocks_t)
  ))
}

# Do the actual simulation
.do_simulate_stdresid_t <- function(t, filtered_copula) {
  u <- .copula_simulate_t(t, filtered_copula)
  garch_uniform2stdresid(model.GARCH, u)
}

# Memoised version for constant 
.do_simulate_stdresid_t_constant <- memoise(
  function(c) .do_simulate_stdresid_t(1, c)
)

do_simulate <- function(name, copula) {
  distribution <- array(NA, dim = c(N_RANDOM, 6, length(T_SEQ)))
  colnames(distribution) <- names(model.GARCH)
  
  for (t_i in seq_along(T_SEQ)) {
    # Select INDEX because we don't want to fill up with all the unsimulated
    # periods
    t <- T_SEQ[t_i]
    
    tic(sprintf('%s: t = %d', name, t))
    distribution[,, t_i] <- garch_path_t(t, copula)
    toc()
  }
  
  filename <- sprintf('%s_%d.RData', name, N_RANDOM)
  file <- file.path('data/derived/distributions', filename)
  
  save(distribution, file = file)
  
  # avoid returning the results of save in case to allow for garbage collection
  return('OK')
}