# Setup ------------------------------------------------------------------

library(devtools)
library(foreach)
library(doParallel)
load_all('australian')

rm(list = ls())

# Gustaf: Create a cluster and pass instead. Remember to stopCluster() after
# i.e.
# cl <- makeCluster()
# registerDoParallel(cl = cl)
#
# It should work!
registerDoParallel(cores = 7)

# Copula Estimation (Full Sample) ----------------------------------------
load('data/derived/garch/model_GARCH_chosen_u.RData')
date <- df.u$Date
df.u <- df.u[, -1]

dynamic_copula <- list(
  norm = CopulaSpecification(
    distribution = CopulaDistribution(nu = Inf, gamma = rep(0, 6)),
    dynamics = CopulaDynamics(alpha = 0.06577153, beta = 0.9145452)
  ),
  
  std = CopulaSpecification(
    distribution = CopulaDistribution(nu = 11.78368353, gamma = rep(0, 6)),
    dynamics = CopulaDynamics(alpha = 0.06881473, beta = 0.91247336)
  ),
  
  ghst = CopulaSpecification(
    distribution = CopulaDistribution(nu = 11.804937402,
                                      gamma = c(
                                        -0.019263030,
                                         0.064795233,
                                        -0.161745377,
                                        -0.148752424,
                                         0.081236747,
                                         0.004100415
                                      )),
    dynamics = CopulaDynamics(alpha = 0.068927164, beta = 0.912195831)
  )
)

norm <- copula_fit(dynamic_copula$norm, df.u, distribution = 'norm', constant = F)
std <- copula_fit(dynamic_copula$std, df.u, distribution = 't', constant = F)
ghst <- copula_fit(dynamic_copula$ghst, df.u, distribution = 'ghst', constant = F)

dynamic_copula_fit <- list(
  norm = norm,
  std = std,
  ghst = ghst
)

save(dynamic_copula_fit, file = 'data/derived/copula/full_dynamic.RData')


# Constant Copula (Full Sample) ------------------------------------------
load('data/derived/garch/model_GARCH_chosen_u.RData')
date <- df.u$Date
df.u <- df.u[, -1]

constant_copula <-list(
  norm = CopulaSpecification(
    distribution = CopulaDistribution(nu = Inf, gamma = rep(0, 6)),
    dynamics = CopulaDynamics(alpha = 0, beta = 0)
  ),
  
  std = CopulaSpecification(
    distribution = CopulaDistribution(nu = 11.80294, gamma = rep(0, 6)),
    dynamics = CopulaDynamics(alpha = 0, beta = 0)
  ),
  
  ghst = CopulaSpecification(
    distribution = CopulaDistribution(nu = 11.80294,
                                      gamma = c(
                                        -0.05473225,
                                         0.07062225,
                                        -0.17046365,
                                        -0.12501372,
                                         0.09492990,
                                         0.02159156
                                      )),
    dynamics = CopulaDynamics(alpha = 0, beta = 0)
  )
)

# XXX Infinite LL for Student's t?
norm <- copula_fit(constant_copula$norm, df.u, distribution = 'norm', constant = T)
std <- copula_fit(constant_copula$std, df.u, distribution = 't', constant = T)
ghst <- copula_fit(constant_copula$ghst, df.u, distribution = 'ghst', constant = T)

constant_copula_fit <- list(
  norm = norm,
  std = std,
  ghst = ghst
)

save(constant_copula_fit, file = 'data/derived/copula/full_constant.RData')

# Copula Estimation (Short Sample) ----------------------------------------

load('data/derived/garch/oos_model_GARCH_chosen_u.RData')
date <- df.u$Date
df.u <- df.u[, -1]

constant_copula <-list(
  norm = CopulaSpecification(
    distribution = CopulaDistribution(nu = Inf, gamma = rep(0, 6)),
    dynamics = CopulaDynamics(alpha = 0, beta = 0)
  ),
  
  std = CopulaSpecification(
    distribution = CopulaDistribution(nu = 8, gamma = rep(0, 6)),
    dynamics = CopulaDynamics(alpha = 0, beta = 0)
  ),
  
  ghst = CopulaSpecification(
    distribution = CopulaDistribution(nu = 8, gamma = rep(0, 6)),
    dynamics = CopulaDynamics(alpha = 0, beta = 0)
  )
)

# XXX Infinite LL for Student's t?
norm <- copula_fit(constant_copula$norm, df.u, distribution = 'norm', constant = T)
std <- copula_fit(constant_copula$std, df.u, distribution = 't', constant = T)
ghst <- copula_fit(constant_copula$ghst, df.u, distribution = 'ghst', constant = T)

constant_copula_fit <- list(
  norm = norm,
  std = std,
  ghst = ghst
)

save(constant_copula_fit, file = 'data/derived/copula/oos_constant.RData')

# Dynamic Copula (Short Sample) ------------------------------------------

load('data/derived/garch/oos_model_GARCH_chosen_u.RData')
date <- df.u$Date
df.u <- df.u[, -1]

dynamic_copula <- list(
  norm = CopulaSpecification(
    distribution = CopulaDistribution(nu = Inf, gamma = rep(0, 6)),
    dynamics = CopulaDynamics(alpha = 0.06, beta = 0.91)
  ),
  
  std = CopulaSpecification(
    distribution = CopulaDistribution(nu = 8, gamma = rep(0, 6)),
    dynamics = CopulaDynamics(alpha = 0.06, beta = 0.91)
  ),
  
  ghst = CopulaSpecification(
    distribution = CopulaDistribution(nu = 8, gamma = rep(0, 6)),
    dynamics = CopulaDynamics(alpha = 0.06, beta = 0.91)
  )
)

norm <- copula_fit(dynamic_copula$norm, df.u, distribution = 'norm', constant = FALSE)
std <- copula_fit(dynamic_copula$std, df.u, distribution = 't', constant = FALSE)
ghst <- copula_fit(dynamic_copula$ghst, df.u, distribution = 'ghst', constant = FALSE)

dynamic_copula_fit <- list(
  norm = norm,
  std = std,
  ghst = ghst
)

save(dynamic_copula_fit, file = 'data/derived/copula/oos_dynamic.RData')

# Filter full datasets ---------------------------------------------------

load('data/derived/garch/oos_model_GARCH_chosen_filtered.RData')
load('data/derived/copula/oos_constant.RData')
load('data/derived/copula/oos_dynamic.RData')

# Compute uniforms from GARCH filtered (this should arguably already be
# done by the fitting thing but whatever, it doesn't take very long)
stdresid <- sapply(filtered, function(f) f@filter$z)
u <- garch_stdresid2uniform(filtered, stdresid)
rm(stdresid)

constant_copula_filtered <- lapply(
  constant_copula_fit,
  function(fit) copula_filter(fit$fit, u)
)

dynamic_copula_filtered <- lapply(
  dynamic_copula_fit,
  function(fit) copula_filter(fit$fit, u)
)

save(constant_copula_filtered, file = 'data/derived/copula/oos_constant_filtered.RData')
save(dynamic_copula_filtered, file = 'data/derived/copula/oos_dynamic_filtered.RData')


# Filter Full Sample with Full Model -------------------------------------

load('data/derived/garch/model_GARCH_chosen_filtered.RData')
load('data/derived/copula/full_dynamic.RData')
load('data/derived/copula/full_constant.RData')

# Compute uniforms from GARCH filtered (this should arguably already be
# done by the fitting thing but whatever, it doesn't take very long)
stdresid <- sapply(filtered, function(f) f@filter$z)
u <- garch_stdresid2uniform(filtered, stdresid)
rm(stdresid)

constant_copula_filtered <- lapply(
  constant_copula_fit,
  function(fit) copula_filter(fit$fit, u)
)

dynamic_copula_filtered <- lapply(
  dynamic_copula_fit,
  function(fit) copula_filter(fit$fit, u)
)

save(constant_copula_filtered, file = 'data/derived/copula/full_constant_filtered.RData')
save(dynamic_copula_filtered, file = 'data/derived/copula/full_dynamic_filtered.RData')

