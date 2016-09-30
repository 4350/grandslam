#' Simulate Return Series from ARMA-GARCH-DCC Copula Model
#' 
#' For each T,
#' 
#' Generate a battery of shocks from copula MV distribution using
#' conditional correlation matrix (for t = 0, use unconditional)
#' 
#' Then, to update correlation matrix for next period:
#' 
#' 1. Standardize shocks
#' 2. Send through DCC model to update correlation matrix for next period
#' 
#' To generate return series:
#' 
#' 1. Perform uniform transformation of (non-standardized) shock according to
#'    copula univariate distributions
#' 2. Compute quantile according to each GARCH model to get GARCH shocks
#' 3. Send through ARMA-GARCH model to get return for this period
#' 
#' Nothing in the copula depends on the GARCH steps; it makes sense to generate
#' the full copula series of shocks first, and then transform these shocks into
#' returns.
#' 
#' The copula simulation is hard to run in parallel, as it involves
#' multivariate recursion.
#' 
#' In order to simulate, the models must be fully specified.
#' 
#' For the copula, I need this:
#' 
#'      list(
#'        dist.params = list(
#'          df = 10,             # NULL if normal
#'          skew = c(3, 5, 6, 3) # NULL if symmetric
#'        ),
#'        alpha = 0.06,
#'        beta = 0.91,
#'        Omega = matrix(c(
#'          1.00, 0.50, 0.50, 0.50,
#'          0.50, 1.00, 0.50, 0.50,
#'          0.50, 0.50, 1.00, 0.50,
#'          0.50, 0.50, 0.50, 1.00
#'        ), nrow = 4, byrow = T)
#'      )
#'  
#'  For each ARMA-GARCH, I need the spec and parameters. This should work
#'  fine with whatever input/output rugarch deals with.

# Setup Libraries ----
library(devtools)
library(rugarch)
library(tictoc)
load_all('wimbledon')
rm(list = ls())

# SETUP ------------------------------------------------------------------
load('data/derived/model_GARCH.RData')
load('data/derived/model_copula_dynamic_ghskt.RData')

kCopulaModel <- model.copula.dynamic.ghskt
kGARCHModels <- model.GARCH
rm(model.GARCH, model.copula.dynamic.ghskt)

T <- 5000
N <- ncol(kCopulaModel$Omega)

# Build copula univariate distributions
sim.c.uv.dists <- function(p.c) {
  p.d <- p.c$dist.params
  
  dist.gauss <- function(i) gauss()
  dist.ght   <- function(i) student.t(nu = p.d$df)
  dist.ghskt <- function(i) student.t(nu = p.d$df, gamma = p.d$skew[[i]])
  
  fn <- dist.gauss
  if (!is.null(p.d$df)) {
    if (!is.null(p.d$skew)) {
      fn <- dist.ghskt
    } else {
      fn <- dist.ght
    }
  }
  
  lapply(seq(ncol(p.c$Omega)), fn)
}

sim.c.rghyp <- function(p.c, Correlation) {
  p.d <- p.c$dist.params
  N <- ncol(p.c$Omega)
  mu <- rep(0, N)
  
  mv.dist <- NULL
  if (is.null(p.d$df)) {
    mv.dist <- gauss(
      mu = mu,
      sigma = Correlation
    )
  }
  else {
    skew <- if (is.null(p.d$skew)) rep(0, N) else p.d$skew
    
    mv.dist <- student.t(
      nu = p.d$df,
      mu = mu,
      gamma = skew,
      sigma = Correlation
    )
  }
  
  ghyp::rghyp(1, mv.dist)
}

sim.c.Q_tp1 <- function(p.c, Q_t, shocks.std_t) {
  (1 - p.c$alpha - p.c$beta) * p.c$Omega +
    p.c$beta * Q_t +
    p.c$alpha * (shocks.std_t %*% t(shocks.std_t))
}

# SIMULATION SETUP -------------------------------------------------------

# Create univariate distribution objects depending on copula parameters
c.uv.dists <- sim.c.uv.dists(kCopulaModel)

# We intialize the correlation to (the correlationalized version of) Omega
Q <- array(dim = c(N, N, T))
Correlation <- NA * Q
Q[,, 1] <- kCopulaModel$Omega
Correlation[,, 1] <- dc.Correlation(array(Q[,, 1], dim = c(N, N, 1)))

shocks <- matrix(ncol = N, nrow = T)
shocks.std <- shocks * NA

# COPULA SIMULATION LOOP -------------------------------------------------

tic()
for (t in seq(T)) {
  # Get shocks from MV distribution for this period
  shocks[t, ] <- sim.c.rghyp(kCopulaModel, Correlation[,, t])
  shocks.std[t, ] <- dc.shocks.std(
    rbind(shocks[t, ]),
    c.uv.dists,
    alpha = kCopulaModel$alpha,
    beta = kCopulaModel$beta
  )
  
  # Prepare Q and Correlation for next period
  if (t < T) {
    Q[,, t + 1] <- sim.c.Q_tp1(
      kCopulaModel,
      Q_t = Q[,, t],
      shocks.std_t = shocks.std[t, ]
    )
    
    Correlation[,, t + 1] <- dc.Correlation(
      array(Q[,, t + 1], dim = c(N, N, 1))
    )
  }
}
toc()
rm(t)

# ARMA-GARCH -------------------------------------------------------------

tic()
garchsim <- sapply(seq(N), function(i) {
  # Compute uniform residuals using the marginal distribution of the copula
  u <- ghyp::pghyp(shocks[, i], c.uv.dists[[i]])
  
  # Compute GARCH shocks from the uniforms using the parametrization of the
  # GARCH models. Note that rugarch uses the GHSKT distribution, however, in
  # a different parametrization to ghyp. However, we want to use ghyp for
  # performance reasons hence the funky function here
  
  # XXX This should use the fit if using a fit!!!
  model <- kGARCHModels[[i]]
  pars <- model@model$pars[, 'Level']
  stdres <- garch.qghyp.rugarch(
    u,
    skew = pars['skew'],
    shape = pars['shape']
  )
  
  # ugarchpath if you have specs
  ugarchsim(
    model,
    n.sim = T,
    m.sim = 1,
    custom.dist = list(name = 'sample', distfit = cbind(stdres))
  )
})
toc()

# SAVE OUTPUT ------------------------------------------------------------

output <- list(
  series = sapply(garchsim, function(p) t(p@simulation$seriesSim)),
  sigma = sapply(garchsim, function(p) t(p@simulation$sigmaSim)),
  resid = sapply(garchsim, function(p) t(p@simulation$residSim))
)

# Hard-Coded Demo Models -------------------------------------------------

# params.copula <- list(
#   dist.params = list(
#     df = 8,
#     skew = c(-0.25, 0.25, 0.10, 0.05, 0.10, -0.05)
#   ),
#   alpha = 0.06,
#   beta = 0.91,
#   Omega = matrix(c(
#     1.00, 0.50, 0.50, 0.50, 0.50, 0.50,
#     0.50, 1.00, 0.50, 0.50, 0.50, 0.50,
#     0.50, 0.50, 1.00, 0.50, 0.50, 0.50,
#     0.50, 0.50, 0.50, 1.00, 0.50, 0.50,
#     0.50, 0.50, 0.50, 0.50, 1.00, 0.50,
#     0.50, 0.50, 0.50, 0.50, 0.50, 1.00
#   ), ncol = 6, byrow = T)
# )
# kGARCHModels <- list(
#   ugarchspec(
#     mean.model = list(armaOrder = c(2, 0)),
#     variance.model = list(
#       model = 'fGARCH',
#       submodel = 'GJRGARCH'
#     ),
#     distribution.model = 'ghst',
#     fixed.pars = list(
#       mu = 0.0011691388,
#       ar1 = 0.75,
#       ar2 = 0.20,
#       omega = 0.0000127421,
#       alpha1 = 0.0851367486,
#       beta1 = 0.8680192797,
#       eta11 = 0.3914889074,
#       skew = -2.5203249430,
#       shape = 12.9541169366
#     )
#   ),
#   ugarchspec(
#     mean.model = list(armaOrder = c(0, 0)),
#     distribution.model = 'ghst',
#     variance.model = list(
#       model = 'fGARCH',
#       submodel = 'GJRGARCH'
#     ),
#     fixed.pars = list(
#       mu = 0.0011691388,
#       omega = 0.0000127421,
#       alpha1 = 0.0851367486,
#       beta1 = 0.8680192797,
#       eta11 = 0.3914889074,
#       skew = -2.5203249430,
#       shape = 12.9541169366
#     )
#   ),
#   ugarchspec(
#     mean.model = list(armaOrder = c(0, 0)),
#     distribution.model = 'ghst',
#     variance.model = list(
#       model = 'fGARCH',
#       submodel = 'GJRGARCH'
#     ),
#     fixed.pars = list(
#       mu = 0.0011691388,
#       omega = 0.0000127421,
#       alpha1 = 0.0851367486,
#       beta1 = 0.8680192797,
#       eta11 = 0.3914889074,
#       skew = -2.5203249430,
#       shape = 12.9541169366
#     )
#   ),
#   ugarchspec(
#     mean.model = list(armaOrder = c(0, 0)),
#     variance.model = list(
#       model = 'fGARCH',
#       submodel = 'GJRGARCH'
#     ),
#     distribution.model = 'ghst',
#     fixed.pars = list(
#       mu = 0.0011691388,
#       omega = 0.0000127421,
#       alpha1 = 0.0851367486,
#       beta1 = 0.8680192797,
#       eta11 = 0.3914889074,
#       skew = -2.5203249430,
#       shape = 12.9541169366
#     )
#   ),
#   ugarchspec(
#     mean.model = list(armaOrder = c(0, 0)),
#     variance.model = list(
#       model = 'fGARCH',
#       submodel = 'GJRGARCH'
#     ),
#     distribution.model = 'ghst',
#     fixed.pars = list(
#       mu = 0.0011691388,
#       omega = 0.0000127421,
#       alpha1 = 0.0851367486,
#       beta1 = 0.8680192797,
#       eta11 = 0.3914889074,
#       skew = -2.5203249430,
#       shape = 12.9541169366
#     )
#   ),
#   ugarchspec(
#     mean.model = list(armaOrder = c(0, 0)),
#     variance.model = list(
#       model = 'fGARCH',
#       submodel = 'GJRGARCH'
#     ),
#     distribution.model = 'ghst',
#     fixed.pars = list(
#       mu = 0.0011691388,
#       omega = 0.0000127421,
#       alpha1 = 0.0851367486,
#       beta1 = 0.8680192797,
#       eta11 = 0.3914889074,
#       skew = -2.5203249430,
#       shape = 12.9541169366
#     )
#   )
# )

