#' Functionality related to GARCH estimation

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
garch.specgen <- function(p, q = 0, r = 1, s = 1, model = 'fGARCH', submodel = 'GJRGARCH',
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

#' rugarch quantile function implemented in ghyp
#'
#' This function translates the rugarch parametrization of the generalized
#' hyperbolic skewed Student's t distribution into the parametrization of
#' the generalized hyperbolic skewed Student's t distribution used by package
#' ghyp, which allows us to estimate using splines which is typically a lot
#' faster for many quantiles
#'
#' @param p univariate quantiles
#' @param shape as specified by rugarch
#' @param skew as specified by rugarch
#'
#' @return vector of quantiles
#' @export
garch.qghyp.rugarch <- function(p, shape, skew) {
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

  ghyp::qghyp(
    p,
    ghyp::student.t(
      mu = mu,
      chi = chi,
      nu = nu,
      sigma = 1,
      gamma = beta
    ),
    method = 'splines'
  )
}

#' Create ugarchspec with fixed.pars from ugarchfit
#' @param fits List of ugarchfit objects
#'
#' @return List of ugarchspec
#' @export
garch.fit2spec <- function(fits) {
  lapply(fits, function(fit) {
    garch.specgen(
      fit@model$modelinc['ar'],
      fit@model$modelinc['ma'],
      fixed.pars = fit@fit$coef,
      vtarget = F # necessary to activate omega
    )
  })
}

#' Function to get values for Standardized Residuals Empirical Density graph
#'
#' @param x fitted model object
#' @return out.df Data frame with histogram bins, limits, normal and sstd lines
#'
#' @export
garch.empirical.density = function(x, ...)
{
  vmodel  = x@model$modeldesc$vmodel
  zseries = as.numeric(residuals(x, standardize=TRUE))
  distribution = x@model$modeldesc$distribution
  idx = x@model$pidx
  pars  = x@fit$ipars[,1]
  skew  = pars[idx["skew",1]]
  shape = pars[idx["shape",1]]
  if(distribution == "ghst") ghlambda = -shape/2 else ghlambda = pars[idx["ghlambda",1]]
  xlim 	= c(min(zseries), max(zseries))
  result 	= hist(x = zseries, col = "grey", border = "white",
                 breaks = "Scott", main = "Empirical Density of Standardized Residuals", xlim = xlim, ylim = c(0,0.6),
                 probability = TRUE, ylab="Probability", cex.main = 0.8,
                 plot = FALSE, ...)
  s = result$mids
  y	= ddist(distribution, s, lambda = ghlambda, skew = skew, shape = shape)
  normline = dnorm(s,0,1)
  histdensity = result$density

  # Output data frame with midpoints of bins (x), normal and sstd line and densities
  out.df <- tidyr::gather(
    data.frame(x = s, Normal = normline, Sstd = y, empirical = histdensity),
    'Distribution',
    'values',
    2:4
  )
}

#' Create nice graphs for diagnostics of ARMA-GARCH fits
#'
#' @param df Data frame - of standardized residuals
#' @param factor String - factor inside the dataframe to consider
#' @param newsimp Data frame - from newsimpact function
#' @param empdens Data frame - from empirical.density function
#'
#' @return null But saves graphs in jpeg format
#'
#' @export
#'
garch.diagplots <- function(df, factor, newsimp, empdens) {
  # Get residuals series
  res <- df[,factor]
  absres <- abs(res)

  newsimpdf <- newsimp[[factor]]
  empdensdf <- empdens[[factor]]

  out.acf <- ggplot2::autoplot(
    acf(res, lag.max = 20, ylab = "", xlab = ""),
    main = "ACF standardized residuals",
    ylab = "",
    xlab = ""
  )+
    theme_Publication()

  out.pacf <- ggplot2::autoplot(
    pacf(res, lag.max = 20, ylab = "", xlab = ""),
    main = "PACF standardized residuals",
    ylim = c(-.1,.1),
    ylab = "",
    xlab = ""
  )+
    theme_Publication()

  out.aacf <- ggplot2::autoplot(
    acf(absres, lag.max = 20),
    main = "ACF standardized absolute residuals",
    ylab = "",
    xlab = ""
  )+
    theme_Publication()

  out.apacf <- ggplot2::autoplot(
    pacf(absres, lag.max = 20, ylab = "", xlab = ""),
    main = "PACF standardized absolute residuals",
    ylim = c(-.1,.1),
    ylab = "",
    xlab = ""
  )+
    theme_Publication()

  out.qq <- ggplot(df, aes_string(sample = factor))+
    stat_qq(distribution = qnorm)+
    scale_y_continuous()+
    ggtitle("QQ plot vs normal distribution")+
    theme_Publication()

  out.ret <- ggplot(df, aes_string(x = 'Date', y = factor))+
    geom_line()+
    xlab("")+
    ylab("")+
    ggtitle("Standardized residuals")+
    coord_cartesian(ylim = c(-5,5))+
    scale_y_continuous()+
    theme_Publication()

  out.newsimpact <- ggplot(newsimpdf, aes(x = x, y = y))+
    geom_line()+
    xlab(expression(epsilon[t - 1]))+
    ylab(expression(sigma[t]^2))+
    coord_cartesian(ylim = c(0,0.025))+
    ggtitle("News impact curve")+
    theme_Publication()


  out.empiricaldensity <- ggplot(empdensdf,
                                 aes(x = x, y = values, group = Distribution)
  ) +
    geom_bar(data = empdensdf[empdensdf$Distribution == c("empirical"), ], stat = 'identity'
    )+
    geom_line(data = empdensdf[empdensdf$Distribution == c("Sstd", "Normal"), ],
              aes(color = Distribution),
              size = 1
    )+
    coord_cartesian(xlim = c(-6, 4), ylim = c(0, .6)) +
    xlab('')+
    ylab('')+
    ggtitle('Empirical density vs distributions')+
    theme_Publication()

  # Grid plots and print to pdf

  g <- arrangeGrob(out.ret, out.qq, out.acf, out.pacf,
                   out.aacf, out.apacf, out.newsimpact, out.empiricaldensity,
                   ncol = 2)
  ggsave(file= paste('output/garch_diagnostics/garch_diagnostics', factor, '.png', sep = ''), g, width = 14.0, height = 21.0, units = 'cm', limitsize = F) #saves g
}
