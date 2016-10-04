#' Functionality for estimating and using a dynamic copula
#'
#' Makes heavy use of parallelization for estimation, epxecting a cluster
#' argument that "knows" about the relevant functions already.

# Series Producing Functions ----
# Functions used to produce the Q series for maximum likelihood estimation
# and simulation.

dc.uv.dists <- function(N, dist.params) {
  lapply(seq(N), function(i) {
    # Gaussian Case
    if (is.null(dist.params$df)) {
      dist <- ghyp::gauss()
    }
    # Symmetric T case
    else if (is.null(dist.params$skew)) {
      dist <- student.t(nu = dist.params$df)
    }
    # Skewed T case
    else {
      dist <- student.t(nu = dist.params$df, gamma = dist.params$skew[i])
    }

    dist
  })
}

#' Invert uniform realizations to the appropriate distributions
#'
#' @param u TxN uniforms
#' @param uv.dists list of N distributions
#' @param cluster parallel cluster to compute on
#'
#' @export
dc.shocks <- function(u, uv.dists, cluster = NULL) {
  invert <- function(i, u, uv.dists) {
    ghyp::qghyp(
      u[, i],
      uv.dists[[i]],

      # Only relevent for asymmetric t; for Gaussian and symmetric t, qghyp
      # automatically calls the appropriate inversion functions
      method = 'splines'
    )
  }

  if (!is.null(cluster)) {
    return(parSapply(cluster, seq(ncol(u)), invert, u = u, uv.dists = uv.dists))
  }

  sapply(seq(ncol(u)), invert, u = u, uv.dists = uv.dists)
}

#' Recursively build the normalized shocks
#'
#' Transform shocks to normalized shocks by dividing through with the
#' diagonal of conditional "correlation" matrix (qdiag).
#'
#' @param shocks
#' @parma uv.dists
#' @param alpha
#' @param beta
#'
#' @export
dc.shocks.std <- function(shocks, uv.dists, alpha, beta) {
  # If innovations are t-distributed, they do not have unit variance; and if
  # skewed, they also don't have expectation zero. We therefore standardize
  # all shocks using the distribution moments. This follows Christoffersen.

  # rbind ensures that this is a matrix even if T = 1
  shocks <- rbind(sapply(seq(ncol(shocks)), function(i) {
    dist <- uv.dists[[i]]
    (shocks[, i] - ghyp::mean(dist)) / sqrt(ghyp::vcov(dist))
  }))

  # Add an initial observation (assumed zero)
  shocks <- rbind(rep(0, ncol(shocks)), shocks)

  # Series that we're building here
  #
  #   shocks.std - shocks that have been divided by the square root of the
  #                conditional variance this period
  #   qdiag      - diagonals of Q_t (conditional variance)
  shocks.std <- shocks * NA
  qdiag <- shocks * NA
  qdiag[1, ] <- 1

  for (t in 2:nrow(shocks)) {
    # Compute diagonal of q for this period according to Christoffersen
    qdiag[t, ] <-
      (1 - alpha - beta) * 1 +
      beta * qdiag[t - 1, ] +
      alpha * (shocks[t - 1, ] / qdiag[t - 1, ]) ^ 2

    # Compute normalized shock this period by dividing with diagonal
    shocks.std[t, ] <- shocks[t, ] / sqrt(qdiag[t, ])
  }

  # Although we computed qdiag already and could in principle save it, we
  # just get it again when we do the real recursion
  # Throw away initial observation
  shocks.std[-1, ]
}

#' Estimate Omega (aka S) as the expectation of the standardized shocks
#'
#' This follows Aielli (2009) and Christoffersen (2012)
#'
#' @param shocks.std copula shocks divided by diagonal of Q_t
#'
#' @export
dc.Omega <- function(shocks.std) {
  Omega <- matrix(0, ncol(shocks.std), ncol(shocks.std))

  # Vectorizing this code does not make it faster but considerably harder
  # to read!
  for (t in 1:nrow(shocks.std)) {
    Omega = Omega + (shocks.std[t, ] %*% t(shocks.std[t, ]))
  }

  Omega / nrow(shocks.std)
}


#' Compute Q series according to cDCC model
#'
#' @param shocks.std standardized/normalized shocks
#' @param Omega long-term expectation of Q
#' @param alpha punch of shocks
#' @param beta autoregressive coefficient
#'
#' @return NxNxT series of Q
#' @export
dc.Q <- function(shocks.std, Omega, alpha, beta) {
  N <- ncol(shocks.std)
  T <- nrow(shocks.std)

  # One extra observation for the first go
  Q <- array(dim = c(N, N, T + 1))
  shocks.std <- rbind(rep(0, N), shocks.std)

  # t0 observation = expectation. Also plausible: diag(1, N)
  Q[,, 1] <- Omega

  for (t in 2:(T + 1)) {
    Q[,, t] <-
      (1 - alpha - beta) * Omega +
      beta * Q[,, t - 1] +
      alpha * (shocks.std[t - 1, ] %*% t(shocks.std[t - 1, ]))
  }

  Q[,, -1]
}

#' Standardize Q to a correlation series
#'
#' @param Q NxNxT matrix
#'
#' @return NxNxT array with unit diagonals
#' @export
dc.Correlation <- function(Q, cluster = NULL) {
  fn <- function(t, Q) {
    Q_t <- Q[,, t]
    inv.sqrt <- solve(sqrt(diag(diag(Q_t))))
    inv.sqrt %*% Q_t %*% inv.sqrt
  }

  seqT <- seq(dim(Q)[3])

  # Run in parallel if a cluster was passed
  if (!is.null(cluster)) {
    return(parSapply(cluster, seqT, fn, Q = Q, simplify = "array"))
  }

  sapply(seqT, fn, Q = Q, simplify = "array")
}

dc.run.model <- function(u, dist.params, alpha, beta, cluster = NULL) {
  # Build univariate distributions as they're used for the construction
  # of our shocks; the MV distribution is built for each t based on the
  # Correlation matrix generated
  uv.dists <- dc.uv.dists(ncol(u), dist.params)
  shocks <- dc.shocks(u, uv.dists, cluster)

  # Get standardized shocks
  shocks.std <- dc.shocks.std(shocks, uv.dists, alpha, beta)

  # Compute Omega using method of moments
  Omega <- dc.Omega(shocks.std)

  # Get correlation time series
  Q <- dc.Q(shocks.std, Omega, alpha, beta)
  Correlation <- dc.Correlation(Q, cluster = cluster)

  list(
    uv.dists = uv.dists,
    shocks = shocks,
    shocks.std = shocks.std,
    Omega = Omega,
    Q = Q,
    Correlation = Correlation
  )
}

# Log Likelihood Functions ----

dc.ll.marginal <- function(shocks, uv.dists, cluster) {
  fn <- function(i, shocks, uv.dists) {
    sum(ghyp::dghyp(shocks[, i], uv.dists[[i]], logvalue = T))
  }
  
  parSapply(
    cluster,
    seq(ncol(shocks)),
    fn,
    
    # Closures on clusters seems to work intermittently; pass all relevant
    # parameters directly for safety (seems to have a slight performance
    # impact).
    shocks = shocks,
    uv.dists = uv.dists
  )
}

dc.ll.joint <- function(shocks, dist.params, Correlation, cluster) {
  # For performance reasons, the skewed t-distribution calls a special density
  # function which skips a lot of checks. This shaves about 500 milliseconds
  # per call relative to calling dghyp.
  skewedTFn <- function(t, shocks, params, Correlation) {
    unname(as.vector(dghsst(
      rbind(shocks[t, ]),
      nu = params$df,
      gamma = params$skew,
      sigma = Correlation[,, t]
    )))
  }

  symmetricTfn <- function(t, shocks, params, Correlation) {
    dist <- ghyp::student.t(
      nu = params$df,
      mu = rep(0, ncol(shocks)),
      sigma = Correlation[,, t]
    )
    ghyp::dghyp(shocks[t, ], dist, logvalue = T)
  }

  gaussianFn <- function(t, shocks, params, Correlation) {
    dist <- ghyp::gauss(
      mu = rep(0, ncol(shocks)),
      sigma = Correlation[,, t]
    )
    ghyp::dghyp(shocks[t, ], dist, logvalue = T)
  }

  fn <- gaussianFn
  if (!is.null(dist.params$df)) {
    if (is.null(dist.params$skew)) {
      fn <- symmetricTfn
    }
    else {
      fn <- skewedTFn
    }
  }
  
  parSapply(
    cluster,
    seq(nrow(shocks)),
    fn,

    # Again: Pass all parameters explicitly seems safer than not doing it
    shocks = shocks,
    params = dist.params,
    Correlation = Correlation
  )
}

dc.ll.total <- function(u, dist.params, alpha, beta, cluster) {
  model <- dc.run.model(
    u = u,
    dist.params = dist.params,
    alpha = alpha,
    beta = beta,
    cluster = cluster
  )
  
  # Compute joint and marginal likelihood
  ll.joint <- sum(dc.ll.joint(
    shocks = model$shocks,
    dist.params = dist.params,
    Correlation = model$Correlation,
    cluster = cluster
  ))
  
  ll.marginal <- sum(dc.ll.marginal(
    shocks = model$shocks,
    uv.dists = model$uv.dists,
    cluster = cluster
  ))
  
  # Sklar's theorem
  ll.joint - ll.marginal
}