#' One-step ahead Q matrix
#'
#' @param p.c copula model
#' @param Q_t Q this period
#' @param shocks.std_t shocks this period
#'
#' @return Q next period
#' @export
sim.c.Q_tp1 <- function(copula, Q_t, shocks.std_t) {
  coef <- copula$params

  (1 - coef$alpha - coef$beta) * copula$Omega +
    coef$beta * Q_t +
    coef$alpha * (shocks.std_t %*% t(shocks.std_t))
}

#' Simulate copula shocks (default 1 period)
#'
#' @param p.c copula parameters
#' @param Correlation Correlation matrix (Sigma)
#' @param n number of shocks (for constant Correlation)
#'
#' @return Nxn shocks
#' @export
sim.c.rghyp <- function(p.c, Correlation, n = 1) {
  p.d <- p.c$params$dist.params
  N <- ncol(p.c$Omega)
  mu <- rep(0, N)

  mv.dist <- NULL
  if (is.null(p.d$df)) {
    mv.dist <- ghyp::gauss(
      mu = mu,
      sigma = Correlation
    )
  }
  else {
    skew <- if (is.null(p.d$skew)) rep(0, N) else p.d$skew

    mv.dist <- ghyp::student.t(
      nu = p.d$df,
      mu = mu,
      gamma = skew,
      sigma = Correlation
    )
  }

  ghyp::rghyp(n, mv.dist)
}

#' Simulate NxT shocks from Copula model
#'
#' @param T number of periods to simulate
#' @param copula parameters describing copula model
#'
#' @return
#' @export
sim.c.copula <- function(T, copula) {
  # Setup
  N <- ncol(copula$Omega)

  # Univariate distributions for shocks
  c.uv.dists <- dc.uv.dists(N, copula$params$dist.params)

  # Copula model results
  shocks <- matrix(ncol = N, nrow = T)
  Q <- array(dim = c(N, N, T))

  shocks.std <- shocks * NA
  Correlation <- NA * Q

  # We intialize the correlation to (the correlationalized version of) Omega
  Q[,, 1] <- copula$Omega
  Correlation[,, 1] <- dc.Correlation(array(Q[,, 1], dim = c(N, N, 1)))

  for (t in seq(T)) {
    # Get shocks from MV distribution for this period
    shocks[t, ] <- sim.c.rghyp(copula, Correlation[,, t])

    shocks.std[t, ] <- dc.shocks.std(
      rbind(shocks[t, ]),
      c.uv.dists,
      alpha = copula$params$alpha,
      beta = copula$params$beta
    )

    # Prepare Q and Correlation for next period
    if (t < T) {
      Q[,, t + 1] <- sim.c.Q_tp1(
        copula,
        Q_t = Q[,, t],
        shocks.std_t = shocks.std[t, ]
      )

      Correlation[,, t + 1] <- dc.Correlation(
        array(Q[,, t + 1], dim = c(N, N, 1))
      )
    }
  }

  # Compute uniform residuals using the marginal distribution of the copula
  sapply(seq(N), function(i) ghyp::pghyp(shocks[, i], c.uv.dists[[i]]))
}

#' Simulate GARCH paths from uniforms (given by copula)
#'
#' @param u NxT uniform residuals
#' @param garch list of N uGARCHfit models
#'
#' @return list of series, sigma and residuals matrices
#' @export
sim.c.GARCH <- function(u, garch) {
  T <- nrow(u)

  sim <- sapply(seq(ncol(u)), function(i) {
    model <- garch[[i]]
    pars <- model@model$pars[, 'Level']

    stdres <- garch.qghyp.rugarch(
      u[, i],
      skew = pars['skew'],
      shape = pars['shape']
    )

    rugarch::ugarchsim(
      model,
      n.sim = T,
      m.sim = 1,
      custom.dist = list(name = 'sample', distfit = cbind(stdres))
    )
  })

  # Collect it all nicely (we don't care about the models or the seed)
  list(
    series = sapply(sim, function(p) t(p@simulation$seriesSim)),
    sigma = sapply(sim, function(p) t(p@simulation$sigmaSim)),
    resid = sapply(sim, function(p) t(p@simulation$residSim))
  )
}

#' Simulate a full series of T observation from GARCH and Copula
#'
#' @param garch List of GARCH models (uGARCHfit objects)
#' @param copula Copula parametrization
#'
#' @return
#' @export
sim.c <- function(T, garch, copula) {
  sim.c.GARCH(sim.c.copula(T, copula), garch)
}

#' Turn mxN copula shocks into mxN standardized residuals
#'
#' @param shocks mxN shocks from copula model
#' @param uv list of univariate distributions
#' @param garch list of GARCH models
#'
#' @return
shocks2stdresid <- function(shocks, uv, garch, cluster = NULL) {
  fn <- function(i, shocks, uv, garch) {
    garch.pars <- garch[[i]]@model$fixed.pars
    copula.dist <- uv[[i]]
    
    # Transform copula shocks into uniforms 
    u <- ghyp::pghyp(shocks[, i], copula.dist)
    
    # Transform uniforms into GARCH stdresid
    wimbledon::garch.qghyp.rugarch(
      u,
      skew = garch.pars['skew'],
      shape = garch.pars['shape']
    )
  }
  
  NN <- seq(ncol(shocks))
  if (is.null(cluster)) {
    return(sapply(NN, fn, garch = garch, uv = uv, shocks = shocks))
  }
  
  parSapply(cluster, NN, fn, garch = garch, uv = uv, shocks = shocks)
}
