#' @export GARCHCopulaModel
GARCHCopulaModel <- setClass('GARCHCopulaModel',
                             slots = c(garch = 'list',
                                       copula = 'CopulaSpecification'))

#' Filter a return series through GARCH Copula Model
#'
#' @param model GARCHCopulaModel
#' @param x Return series
#' @param X Upsilon regressors
#'
#' @return result of copula_filter
#' @export
garch_copula_filter <- function(model, x, X = NULL) {
  # GARCH filter to uniforms
  u <- foreach(i = seq_along(model@garch), .combine = 'cbind') %do% {
    stdresid <- ugarchfilter(model@garch[[i]], x[, i])@filter$z

    pars <- model@garch[[i]]@model$pars[, 'Level']

    rugarch:::psghst(stdresid,
      shape = pars['shape'],
      skew = pars['skew']
    )
  }

  copula_filter(spec, u, X)
}
