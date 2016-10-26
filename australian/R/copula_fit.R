#' Fit copula to uniform residuals.
#'
#' Note: This takes a long time.
#'
#' @param spec CopulaSpecification to fit
#' @param u TxN uniform residuals
#' @param NxkxT array of external regressors
#' @param constant Optimize constant copula (alpha = beta = 0)
#'
#' @return Fitted CopulaSpecification
#' @export
copula_fit <- function(spec, u, X = NULL, distribution = 'norm',
                       constant = TRUE, upsilon = FALSE) {
  # Figure out what parameters we're optimizing over
  dynamics <- spec@dynamics
  N <- ncol(u)

  if (constant && !all(c(dynamics@alpha, dynamics@beta) == 0)) {
    stop('Constant copula requires alpha = beta = 0')
  }

  # Can't do Upsilon optimization with constant copula or if no X
  stopifnot(!(constant && upsilon))
  stopifnot(!(upsilon && is.null(X)))

  # Gets appropriate theta
  pars <- copula_fit_pars(spec, distribution, constant, upsilon)

  # Note: Constant normal copula has no parameters; but run a single optim
  # anyway just to get optimized value. constrOptim doesn't work with empty
  # ui/ci

  if (length(pars$pars) == 0) {
    optimized <- stats::optim(pars$pars, copula_optimize,
                              spec = spec, u = u, X = X,
                              distribution = distribution,
                              constant = constant,
                              upsilon = upsilon)
  }
  else if (length(pars$pars) == 1) {
    stopifnot(length(pars$ci) == 2)
    stopifnot(distribution == 't')
    
    optimized <- stats::optim(
      pars$pars,
      copula_optimize,
      
      spec = spec,
      u = u,
      X = X,
      distribution = distribution,
      constant = constant,
      upsilon = upsilon,
      
      method = 'Brent',
      lower = pars$ci[1],
      upper = -pars$ci[2]
    )
  }
  else {
    optimized <- stats::constrOptim(
      pars$pars,
      copula_optimize,
      grad = NULL,

      spec = spec,
      u = u,
      X = X,
      distribution = distribution,
      constant = constant,
      upsilon = upsilon,

      ui = pars$ui,
      ci = pars$ci,
      
      control = list(
        trace = 6,
        maxit = 1000
      )
    )
  }

  # Rebuild and filter the final specification (to get Omega)
  fit <- copula_fit_build_spec(optimized$par, spec, distribution, constant, upsilon)
  filtered <- copula_filter(fit, u, X)

  list(
    fit = filtered$spec,
    ll = filtered$ll,
    optimized = optimized
  )
}

copula_optimize <- function(pars, spec, u, X, distribution, constant, upsilon) {
  spec <- copula_fit_build_spec(pars, spec, distribution, constant, upsilon)
  -copula_filter(spec, u, X)$ll
}

#' Get pars and constraints for a copula model
#'
#' @param N how many factors?
#' @param distribution norm, t or ghst
#' @param constant constant copula?
#' @param upsilon with upsilon?
#'
#' @return list with theta, ui and ci
#' @export
copula_fit_pars <- function(spec, distribution = 'norm', constant = T, upsilon = F) {
  N <- length(spec@distribution@gamma)

  if (upsilon && !all(spec@dynamics@theta == spec@dynamics@theta[1])) {
    stop('Only support for single theta')
  }

  # Distribution boundaries
  INITIAL_NU <- spec@distribution@nu
  INITIAL_GAMMA <- spec@distribution@gamma
  UI_NU <- rbind(1, -1)
  CI_NU <- rbind(6, -20)
  GAMMA_B <- 0.25

  # Initial and restrictions on alpha and beta
  INITIAL_ALPHABETA <- c(spec@dynamics@alpha, spec@dynamics@beta)
  UI_ALPHABETA <- rbind(diag(1, 2), c(-1, -1))
  CI_ALPHABETA <- rbind(0, 0, -0.9999)

  # Initial value and restrictions on phi (with Upsilon)
  INITIAL_PHI <- spec@dynamics@phi
  UI_PHI <- rbind(1, -1)
  CI_PHI <- rbind(0, -1)

  # Only support a single theta
  INITIAL_THETA <- spec@dynamics@theta[1,1]

  pars <- list()
  ui <- list()
  ci <- list()

  # Distribution parameters to optimize over
  if (distribution == 't') {
    pars$distribution <- c(INITIAL_NU)
    names(pars$distribution) <- c('nu')

    ui$distribution <- UI_NU
    ci$distribution <- CI_NU
  }
  else if (distribution == 'ghst') {
    pars$distribution <- c(INITIAL_NU, INITIAL_GAMMA)
    names(pars$distribution) <- c('nu', paste0('gamma', 1:N))

    # Add gamma restrictions
    ui$distribution <- magic::adiag(UI_NU,
                                    rbind(diag(1, N), diag(-1, N)))
    ci$distribution <- rbind(CI_NU,
                             cbind(rep(-GAMMA_B, N)), # lower
                             cbind(rep(-GAMMA_B, N))) # uppper
  }

  if (!constant) {
    pars$alphabeta <- INITIAL_ALPHABETA
    names(pars$alphabeta) <- c('alpha', 'beta')

    ui$alphabeta <- UI_ALPHABETA
    ci$alphabeta <- CI_ALPHABETA
  }

  if (upsilon) {
    pars$upsilon <- c(INITIAL_PHI, INITIAL_THETA)
    names(pars$upsilon) <- c('phi', 'theta')

    # theta is unrestricted
    ui$upsilon <- cbind(UI_PHI, 0)
    ci$upsilon <- CI_PHI
  }

  if (length(ui)) {
    ui <- do.call(magic::adiag, ui)
  }
  else {
    ui <- NULL
  }

  list(
    pars = unlist(pars),
    ui = ui,
    ci = unlist(ci)
  )
}

copula_fit_build_spec <- function(pars, spec, distribution, constant, upsilon) {
  N <- length(spec@distribution@gamma)

  if (distribution == 't') {
    spec@distribution@nu <- pars[1]
    pars <- pars[-1]
  }
  else if (distribution == 'ghst') {
    spec@distribution@nu <- pars[1]
    spec@distribution@gamma <- pars[2:(N + 1)]
    pars <- pars[-1:-(N + 1)]
  }

  if (!constant) {
    spec@dynamics@alpha <- pars[1]
    spec@dynamics@beta <- pars[2]
    pars <- pars[-1:-2]
  }

  if (upsilon) {
    spec@dynamics@phi <- pars[1]
    spec@dynamics@theta <- cbind(rep(pars[2], N))
    pars <- pars[-1:-2]
  }

  # We should have cleared all pars by now
  stopifnot(length(pars) == 0)

  spec
}
