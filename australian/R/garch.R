#' @importFrom rugarch ugarchfit
NULL

#' Fit GARCH models to data
#'
#' @param garch List of ugarchspec
#' @param data TxN data points
#'
#' @return
#' @export
garch_fit <- function(garch, data) {
  foreach(i = seq_along(garch)) %dopar% {
    ugarchfit(garch[[i]], data[, i], solver = 'hybrid')
  }
}

#' Turns GARCH fits to specs with fixed.pars
#'
#' @param fits List of GARCH fits
#'
#' @return List of GARCH specs
#' @export
garch_fit2spec <- function(fits) {
  lapply(fits, function(fit) {
    garch_specgen(
      fit@model$modelinc['ar'],
      fit@model$modelinc['ma'],
      fixed.pars = fit@fit$coef,
      vtarget = F # necessary to activate omega
    )
  })
}

#' Create GARCH specifications
#'
#' Makes building GARCH models easier as standard options
#' apart from p, q order is preset. Standard options include
#' variance targeting, GJRGACH(1,1) specification and the GHST
#' distribution for innovations
#'
#' @param p Scalar - AR order
#' @param q Scalar - MA order
#' @param r Scalar - GARCH alpha order
#' @param s Scalar - GARCH beta order
#' @param vtarget Logical - setting long-term variance using unconditional mean
#'
#' @return uGARCHspec object
#' @export
garch_specgen <- function(p, q = 0, r = 1, s = 1, model = 'fGARCH', submodel = 'GJRGARCH',
                          vtarget = T, dist = 'ghst', fixed.pars = list()) {
  spec <- rugarch::ugarchspec(
    mean.model = list(
      armaOrder = c(p,q)
    ),
    distribution.model = dist,
    variance.model = list(
      model = model,
      submodel = submodel,
      variance.targeting = vtarget
    ),
    fixed.pars = fixed.pars
  )

  spec
}

#' Turn uniform to stdresid for a list of GARCH models
#'
#' @param garch N list of GARCH models
#' @param u TxN uniforms
#'
#' @export
garch_uniform2stdresid <- function(garch, u) {
  foreach (i = seq_along(garch), .combine = 'cbind') %dopar% {
    .garch_qghyp(garch[[i]], u[, i])
  }
}

#' Turn standardized residuals to uniforms. Only works with ghst GARCH
#'
#' @param garch N list of GARCH models
#' @param stdresid TxN matrix of stdresid
#'
#' @export
garch_stdresid2uniform <- function(garch, stdresid) {
  foreach(i = seq_along(garch), .combine= 'cbind') %do% {
    pars <- garch[[i]]@model$pars[, 'Level']

    rugarch:::psghst(stdresid[, i], shape = pars['shape'], skew = pars['skew'])
  }
}

.garch_qghyp <- function(garch_i, p) {
  pars <- garch_i@model$pars[, 'Level']

  shape <- pars['shape']
  skew <- pars['skew']

  # The rugarch parametrization is *kinda* the location and scale invariant
  # parametrization mentioned in the ghyp package documentation. The below
  # code is copy-pasted from rugarch source code (rugarch-distributions.R)
  # which uses the SkewHyperbolic library, and then fed to the ghyp
  # quantile function.
  #
  # chi is delta ^ 2
  nu <- shape
  chi <- 1 / ( ((2 * skew^2)/((nu-2)*(nu-2)*(nu-4))) + (1/(nu-2)) )
  beta <- skew / sqrt(chi)
  mu <- -(beta * chi / (nu - 2))

  ghyp::qghyp(p, ghyp::student.t(
    mu = mu,
    chi = chi,
    nu = nu,
    sigma = 1,
    gamma = beta
  ), method = 'splines')
}
