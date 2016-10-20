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

# QQ plots ----------------------------------------------------------------

#' Gets the data for the qq plot for all factors best fits
#' Based on rugarch but uses faster code for ghst quantile
#' 
do_qq_data <- function(model.GARCH) {
  
  # Get list of qq data for all fits
  qq_data <- lapply(model.GARCH, function(fit) garch.qq(fit))
  
  # Give levels for facet
  qq_data <- bind_rows(qq_data, .id = 'factor')
  qq_data$order <- factor(qq_data$factor, levels = names(model.GARCH))
  return(qq_data)
}

do_qq_plot <- function(qq_data) {
  # Run ggplot on the qq_data
  ggplot(qq_data, aes(x = theo_x, y = sample_y)) +
    geom_point(size = 1) +
    geom_abline(linetype = 2, intercept = 0, slope = 1) +
    facet_wrap(~ order, nrow = 2, ncol = 3)+
    theme_Publication() +
    ylab('Sample quantile')+
    xlab('Theoretical quantile')
  
}

#' Function to get values for QQ plot
#'
#' @param x fitted model object
#' @return out.df Data frame with quantiles for sample and theoretical
#'
#' @export
garch.qq = function(x, ...)
{
  vmodel  = x@model$modeldesc$vmodel
  zseries = as.numeric(residuals(x, standardize=TRUE))
  distribution = x@model$modeldesc$distribution
  idx = x@model$pidx
  pars  = x@fit$ipars[,1]
  skew  = pars[idx["skew",1]]
  shape = pars[idx["shape",1]]
  if(distribution == "ghst") ghlambda = -shape/2 else ghlambda = pars[idx["ghlambda",1]]
  out.df <- .qqDist_re(y = zseries, dist = distribution, lambda = ghlambda, skew = skew, shape = shape)
}


#' Quantile-Quantile supporting functions from rugarch, reworked to return x, y
#' and use special ghst qghyp function for speed
.qqDist_re = function (y, dist = "norm", ylim = NULL, main = paste(dist, "- QQ Plot"),
                    xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", doplot = FALSE,
                    datax = FALSE, cex.main = 0.8, shape = NULL, skew = NULL, ...)
{	
  y = as.vector(y)
  if (has.na <- any(ina <- is.na(y))) {
    yN = y
    y = y[!ina]
  }
  if (0 == (n <- length(y))) stop("y is empty or has only NAs")
  if (dist == "ghst") {
    x = garch.qghyp.rugarch(p = ppoints(n), shape = shape, skew = skew)[order(order(y))]
  } else {
    x = qdist(distribution = dist, p = ppoints(n), ...)[order(order(y))]
  }
  
  
  if (has.na) {
    y = x
    x = yN
    x[!ina] = y
    y = yN
  }
  
  data.frame(theo_x = x, sample_y = y)
}



# End table GARCH diagnostics ---------------------------------------------

do_garch_end_table <- function(object) {
  # Count obs and get LLH and BIC
  obs = length(object@fit$residuals)
  LLH = object@fit$LLH
  BIC = tryCatch(infocriteria(object)['Bayes',], error = function(err) NA)
  # Get standardized residuals and dof
  stdresid = object@fit$residuals/object@fit$sigma
  modelinc = object@model$modelinc
  df = sum(modelinc[2:3])
  # Calculate LB tests
  LB_5 = Weighted.Box.test(stdresid, lag = 5, type = "Ljung-Box", fitdf = df)$p.value
  LB_10 = Weighted.Box.test(stdresid, lag = 10, type = "Ljung-Box", fitdf = df)$p.value
  # Calculate ARCH LM tests
  gdf = sum(modelinc[8:9])
  LM_5 = Weighted.LM.test(as.vector(residuals(object)), sigma(object)^2, lag = 5, type = "correlation", fitdf = gdf)$p.value
  LM_10 = Weighted.LM.test(as.vector(residuals(object)), sigma(object)^2, lag = 10, type = "correlation", fitdf = gdf)$p.value
  # Calculate sign bias tests
  SB_negative = .signbiasTest_re(object)[1]
  SB_positive = .signbiasTest_re(object)[2]
  # Put together out vector
  out.v <- c(
    obs = obs,
    LLH = LLH,
    BIC = BIC,
    LB_5 = LB_5,
    LB_10 = LB_10,
    LM_5 = LM_5,
    LM_10 = LM_10,
    SB_negative = SB_negative,
    SB_positive = SB_positive
  ) %>% round(digits = 4)
  return(out.v)
}

#' Adapted from rugarch to only output p-values for positive and negative
.signbiasTest_re = function(object)
{
  if(is(object, "uGARCHfilter")) z = object@filter$z else z = z = object@fit$z
  res = as.numeric(residuals(object))
  z2 = z^2
  n = length(z)
  zminus = as.integer(res<0)
  zplus = 1-zminus
  czminus = zminus*res
  czplus = zplus*res
  cz = cbind(rep(1, n), zminus, czminus, czplus)
  cz = cz[1:(n-1),]
  z2 = matrix(z2[2:n],ncol = 1)
  cm = data.frame(y = z2, const = cz[,1], zminus = cz[,2], czminus = cz[,3], czplus = cz[,4])
  fitA = lm(y~const+zminus+czminus+czplus-1, data = cm)
  resA = residuals(fitA)
  sgmatrix = matrix(ncol = 3, nrow = 4)
  rownames(sgmatrix) = c("Sign Bias","Negative Sign Bias","Positive Sign Bias","Joint Effect")
  colnames(sgmatrix) = c("t-value","prob","sig")
  sgmatrix[1:3,1:2] = abs(summary(fitA)$coef[2:4,3:4])
  jeffect = rugarch:::.linear.hypothesis(fitA, c("zminus = 0", "czminus = 0","czplus  =  0"), test  =  "Chisq")
  sgmatrix[4,1] = jeffect[5]$Chisq[2]
  sgmatrix[4,2] = jeffect[6][2,]
  sgmatrix = as.data.frame(sgmatrix)
  sgmatrix[,3] = rugarch:::.stars(sgmatrix[,2])
  return(sgmatrix[2:3,2])
}

#  ------------------------------------------------------------------------



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
