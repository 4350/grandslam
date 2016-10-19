#' @importFrom foreach foreach %do% %dopar%
#' @importFrom methods new
NULL

setClassUnion("matrixOrNULL", members = c("matrix", "NULL"))

#' @export CopulaDistribution
CopulaDistribution <- setClass('CopulaDistribution',
                               slots = c(nu = 'numeric',
                                         gamma = 'numeric'),
                               prototype = list(nu = Inf))

#' @export CopulaDynamics
CopulaDynamics <- setClass('CopulaDynamics',
                           slots = c(alpha = 'numeric',
                                     beta = 'numeric',
                                     phi = 'numeric',
                                     theta = 'matrixOrNULL',
                                     Omega = 'matrixOrNULL'),
                           prototype = list(alpha = 0, beta = 0, phi = 0,
                                            theta = NULL, Omega = NULL))

#' @export CopulaSpecification
CopulaSpecification <- setClass('CopulaSpecification',
                                slots = c(distribution = 'CopulaDistribution',
                                          dynamics = 'CopulaDynamics'))

#' @export
copula_filter <- function(spec, u, X = NULL) {
  shocks <- .copula_shocks(spec, u)

  # Standardize and conditionalize shocks as appropriate for DCC model. The
  # naming is kind of unfortunate here; legacy from before.
  shocks_std <- .copula_shocks_standardize(spec, shocks)
  shocks_std <- .copula_shocks_conditionalize(spec, shocks_std)

  if (is.null(X)) {
    # phi must be 0 if you did not supply X
    stopifnot(spec@dynamics@phi == 0)

    # Create a dummy series of Upsilon, so that we don't have to check for
    # it in every function. Note: We don't create an NxN array because no
    # function should care about Upsilon if phi == 0, so just save space and
    # do 1x1
    Upsilon <- array(0, c(1, 1, nrow(u)))
  }
  else {
    Upsilon <- copula_Upsilon(spec@dynamics@theta, X)
  }

  if (is.null(spec@dynamics@Omega)) {
    spec@dynamics@Omega <- .copula_Omega(spec, shocks_std, Upsilon)
  }

  Correlation <- .copula_Correlation(.copula_Q(spec, shocks_std, Upsilon))

  scores <- .copula_scores(spec, shocks, Correlation)

  list(
    spec = spec,
    shocks = shocks,
    Correlation = Correlation,
    scores = scores,
    ll = sum(scores)
  )
}

#' Simulate uniform residuals from the copula
#'
#' @param spec CopulaSpecification
#' @param n.sim Simulation horizon
#' @param m.sim Number of simulation
#' @param Q_T The final Q_t your dataset (default: Omega)
#' @param shocks_T Final shocks in dataset
#' @param X Nxkxn.sim exogenous regressors
#'
#' @return List of m.sim n.simxN uniform residuals
#' @export
copula_simulate <- function(spec, n.sim, m.sim, Q_T = NULL, shocks_T = NULL, X = NULL) {
  if (is.null(Q_T)) {
    stopifnot(!is.null(spec@dynamics@Omega))

    if (!copula_is_constant(spec)) {
      warning("You are simulating from a dynamic copula but did not provide ",
              "your final Q_t matrix. This is probably a mistake. Omega will ",
              "be assumed")
    }

    Q_T <- spec@dynamics@Omega
  }

  if (is.null(shocks_T)) {
    if (!copula_is_constant(spec)) {
      warning("You are simulating from a dynamic copula but did not provide ",
              "your final shocks. This is probably a mistake. We will assume ",
              "zero shocks at time T")
    }

    shocks_T <- rbind(rep(0, ncol(Q_T)))
  }

  if (is.null(X)) {
    stopifnot(spec@dynamics@phi == 0)
    Upsilon <- array(0, c(1, 1, n.sim))
  }
  else {
    Upsilon <- copula_Upsilon(spec@dynamics@theta, X)
  }

  uv_distributions <- .copula_uv_distributions(spec)

  foreach(i = 1:m.sim) %do% {
    shocks <- .copula_simulate_shocks(spec, n.sim, Q_T, shocks_T, Upsilon)

    # Compute the uniform residuals using the marginal distributions of
    # the copula
    sapply(seq_along(uv_distributions),
           function(i) ghyp::pghyp(shocks[, i], uv_distributions[[i]]))
  }
}

#' @export
copula_is_constant <- function(spec) {
  all(c(spec@dynamics@alpha, spec@dynamics@beta) == 0)
}

#' Build Upsilon array from Nxk coefficients and independent variables
#'
#' It's impossible for the copula to know how many of your theta are "free"
#' or how many of your X are identical; it just takes the Upsilon series as
#' given. You can construct it yourself, or use this function to build an
#' appropriate form following the appendix in Christoffersen et al. (2012).
#'
#' @param theta Nxk matrix of regressors
#' @param X Nxkxt array of independent data
#'
#' @return Nxkxt Upsilon array (correlation matrix)
#' @export
copula_Upsilon <- function(theta, X) {
  stopifnot(ncol(theta) == ncol(X))
  stopifnot(nrow(theta) == nrow(X))

  N <- nrow(X)
  k <- ncol(X)
  T <- dim(X)[3]

  Upsilon <- array(NA, c(N, N, T))

  for (t in seq(T)) {
    A_t <- cbind(diag(1, N), theta * X[,,t])

    # Compute each row's root mean square; then dividing by this vector
    # will actually work as if we did it row-by-row; R recycles the vector
    # as necessary and works column-by-column
    row_RMS <- sqrt(rowSums(A_t ^ 2))
    A_bar <- A_t / row_RMS

    Upsilon[,, t] <- A_bar %*% t(A_bar)
  }

  Upsilon
}

.copula_rghyp <- function(n, spec, Correlation) {
  ghyp::rghyp(n, .copula_mv_distribution(spec, Correlation))
}


# Distributions ----------------------------------------------------------

.copula_uv_distributions <- function(spec) {
  lapply(spec@distribution@gamma, function(gamma_i) {
    if (is.infinite(spec@distribution@nu)) {
      return(ghyp::gauss())
    }

    ghyp::student.t(nu = spec@distribution@nu, gamma = gamma_i)
  })
}

.copula_mv_distribution <- function(spec, Correlation) {
  mu <- rep(0, ncol(Correlation))

  if (is.infinite(spec@distribution@nu)) {
    return(ghyp::gauss(mu = mu, sigma = Correlation))
  }

  ghyp::student.t(mu = mu, nu = spec@distribution@nu,
                  gamma = spec@distribution@gamma, sigma = Correlation)
}

# Filtering --------------------------------------------------------------

#' Compute shocks in the copula model
#'
#' @param spec A CopulaSpecification
#' @param u TxN matrix of uniforms
#'
#' @return shocks TxN matrix of distributed shocks
.copula_shocks <- function(spec, u) {
  uv_distributions <- .copula_uv_distributions(spec)

  shocks <- foreach(i = 1:ncol(u)) %dopar% {
    ghyp::qghyp(u[, i], uv_distributions[[i]], method = 'splines')
  }

  matrix(unlist(shocks), ncol = ncol(u))
}

#' If shocks are not Gaussian, they don't have expectation zero and unit
#' variance, so we subtract the distribution mean and divide by distribution
#' standard deviation
#'
#' @param spec CopulaSpecification
#' @param shocks TxN shocks
#'
#' @return TxN standardized shocks
.copula_shocks_standardize <- function(spec, shocks) {
  uv_distributions <- .copula_uv_distributions(spec)
  shocks <- rbind(shocks)

  standardize <- function(dist, i)
    (shocks[, i] - ghyp::mean(dist)) / sqrt(ghyp::vcov(dist))

  rbind(mapply(standardize, uv_distributions, 1:ncol(shocks)))
}

#' Conditionalize shocks by dividing them with the conditional standard
#' deviation in th DCC model. See Christoffersen for details and whatever
#' paper he refers to for justification.
#'
#' This function is only called as part of the recursion to build Omega.
#' It should not be called if you already have Q_t; just divide by
#' sqrt(diag(Q_t)) (see simulation code)
#'
#' Note: This is different from standardizing which should already be done!
#' I.e. you give me standard UNCONDITIONAL shocks and I give you CONDITIONAL
#' shocks.
#'
#' @param spec CopulaSpecification
#' @param shocks TxN shocks (standardized)
#'
#' @return TxN matrix of conditional shocks
.copula_shocks_conditionalize <- function(spec, shocks) {
  T <- nrow(shocks)

  # Step 2: Following Christoffersen's first paper, we also scale the shocks
  # by the "conditional variance" for this period. Supposedly improves the
  # consistency of estimates.
  alpha <- spec@dynamics@alpha
  beta <- spec@dynamics@beta

  # Add an initial observation, assumed zero and set initial diagonal of Q = 1
  shocks <- rbind(0, shocks)
  shocks_std <- shocks * NA
  qdiag <- shocks * NA
  qdiag[1, ] <- 1

  for (t in 2:(T + 1)) {
    qdiag[t, ] <-
      (1 - alpha - beta) +
      beta * qdiag[t - 1, ] +
      alpha * (shocks[t - 1, ] / qdiag[t - 1, ]) ^ 2

    shocks_std[t, ] <- shocks[t, ] / sqrt(qdiag[t, ])
  }

  # Drop the initial observation
  shocks_std[-1, ]
}

.copula_Omega <- function(spec, shocks_std, Upsilon, use_cor = F) {
  # Christoffersen does write out the second formula, but this might be a
  # sloppy way of meaning the sample correlation; difference should be minimal
  # since shocks should be standardized, however, the sample might be "off"
  if (use_cor) {
    correlation <- stats::cor(shocks_std)
  }
  else {
    correlation <- (t(shocks_std) %*% shocks_std) / nrow(shocks_std)
  }

  if (spec@dynamics@phi == 0) {
    return(correlation)
  }

  # We have to deal with Upsilon; compute element-wise average and weight
  # together to Omega
  mean_Upsilon <- apply(Upsilon, c(1, 2), mean)
  phi <- spec@dynamics@phi

  (correlation - phi * mean_Upsilon) / (1 - phi)
}

.copula_Q_tp1 <- function(spec, Q_t, shocks_std_t, Upsilon_t) {
  alpha <- spec@dynamics@alpha
  beta <- spec@dynamics@beta
  Omega <- spec@dynamics@Omega
  phi <- spec@dynamics@phi

  stopifnot(!is.null(Omega))

  (1 - alpha - beta) * ((1 - phi) * Omega + phi * Upsilon_t) +
    beta * Q_t +
    alpha * shocks_std_t %*% t(shocks_std_t)
}

.copula_Q <- function(spec, shocks_std, Upsilon) {
  N <- ncol(shocks_std)
  T <- nrow(shocks_std)

  # One extra observation at the start
  Q <- array(dim = c(N, N, T + 1))
  shocks_std <- rbind(0, shocks_std)

  # t0 observation = unconditional; diag(1, N) is also plausible
  Q[,, 1] <- spec@dynamics@Omega

  for (t in 2:(T + 1)) {
    # Upsilon_t is contemporaneous with Q_t, however, we don't add an
    # additional row to it so we still pick the element from t - 1
    Upsilon_t <- Upsilon[,, t - 1]
    Q[,, t] <- .copula_Q_tp1(spec, Q[,, t - 1], shocks_std[t - 1, ], Upsilon_t)
  }

  Q[,, -1]
}

.copula_Correlation <- function(Q) {
  fn <- function(t) {
    Q_t <- Q[,, t]
    inv_sqrt <- solve(sqrt(diag(diag(Q_t))))
    inv_sqrt %*% Q_t %*% inv_sqrt
  }

  N <- dim(Q)[1]
  T <- dim(Q)[3]

  correlation_list <- lapply(seq(T), fn)

  array(unlist(correlation_list), dim = c(N, N, T))
}

# Simulation -------------------------------------------------------------

.copula_simulate_shocks <- function(spec, n.sim, Q_T, shocks_T, Upsilon) {
  N <- ncol(Q_T)

  # If the copula is constant, there is no need to update the dynamics.
  # We can just simulate n.sim shocks directly
  if (copula_is_constant(spec)) {
    Correlation <- .copula_Correlation(array(Q_T, dim = c(N, N, 1)))
    shocks <- .copula_rghyp(n.sim, spec, Correlation[,, 1])
    return(shocks)
  }

  # We will simulate n.sim shocks but we need to update the dynamics using
  # the end of the current series, so add one row for that
  shocks <- matrix(ncol = N, nrow = n.sim + 1)
  shocks_std <- shocks
  Q <- array(dim = c(N, N, n.sim + 1))
  Correlation <- Q

  # Initiate dynamics with the end of the series (a model of higher order
  # might support Q_T being an array like preresiduals for rugarch)
  shocks[1, ] <- shocks_T
  shocks_std[1, ] <- .copula_shocks_standardize(spec, shocks[1, ]) / sqrt(diag(Q_T))
  Q[,, 1] <- Q_T
  Correlation[,, 1] <- .copula_Correlation(array(Q[,, 1], dim = c(N, N, 1)))

  for (t in 2:(n.sim + 1)) {
    # Update Q as before; Note that we pick Upsilon[,, t - 1] as in the filter
    # because we didn't add an extra row to it.
    Q[,, t] <- .copula_Q_tp1(
      spec,
      Q[,, t - 1],
      shocks_std[t - 1, ],
      Upsilon[,, t - 1]
    )
    Correlation[,, t] <- .copula_Correlation(array(Q[,, t], c(N, N, 1)))

    shocks[t, ] <- .copula_rghyp(1, spec, Correlation[,, t])
    shocks_std[t, ] <- .copula_shocks_standardize(spec, shocks[t, ]) / sqrt(diag(Q[,, t]))
  }

  rbind(shocks[-1, ])
}

# Log Likelihoods --------------------------------------------------------

.copula_scores <- function(spec, shocks, Correlation) {
  joint <- .copula_ll_joint(spec, shocks, Correlation)
  marginal <- .copula_ll_marginal(spec, shocks)

  joint - marginal
}

.copula_ll_marginal <- function(spec, shocks) {
  uv_distributions <- .copula_uv_distributions(spec)

  fn <- function(dist, i) ghyp::dghyp(shocks[, i], dist, logvalue = TRUE)

  ll <- foreach(i = 1:ncol(shocks), .combine = 'cbind') %do% {
    ghyp::dghyp(shocks[, i], uv_distributions[[i]], logvalue = TRUE)
  }

  cbind(rowSums(ll))
}

.copula_ll_joint <- function(spec, shocks, Correlation) {
  fn <- function(t) {
    mvtnorm::dmvnorm(shocks[t, ], sigma = Correlation[,, t], log = TRUE)
  }

  # Use our special density function that cuts of a lot of the fat
  if (is.finite(spec@distribution@nu)) {
    fn <- function(t) {
      temp <- .dghst(rbind(shocks[t, ]),
              nu = spec@distribution@nu,
              gamma = spec@distribution@gamma,
              sigma = Correlation[,, t])

      unname(as.vector(temp))
    }
  }

  # If you ever feel the need to verify our own density function, uncomment
  # the below lines to use ghyp package. It takes CONSIDERABLY longer.
  # fn <- function(t) {
  #   mv_distribution <- .copula_mv_distribution(spec, Correlation[,, t])
  #   ghyp::dghyp(shocks[t, ], mv_distribution, logvalue = TRUE)
  # }

  # Running this in parallel actually seems to have little benefit
  ll <- lapply(seq(nrow(shocks)), fn)
  cbind(unlist(ll))
}
