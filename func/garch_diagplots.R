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
  )
  
  out.pacf <- ggplot2::autoplot(
    pacf(res, lag.max = 20, ylab = "", xlab = ""),
    main = "PACF standardized residuals",
    ylim = c(-.1,.1),
    ylab = "", 
    xlab = ""
  )
  
  out.aacf <- ggplot2::autoplot(
    acf(absres, lag.max = 20),
    main = "ACF standardized absolute residuals",
    ylab = "", 
    xlab = ""
  )
  
  out.apacf <- ggplot2::autoplot(
    pacf(absres, lag.max = 20, ylab = "", xlab = ""),
    main = "PACF standardized absolute residuals",
    ylim = c(-.1,.1),
    ylab = "", 
    xlab = ""
  )
  
  out.qq <- ggplot(df, aes_string(sample = factor))+ 
    stat_qq(distribution = qnorm)+
    scale_y_continuous()+
    ggtitle("QQ plot vs normal distribution")
  
  out.ret <- ggplot(df, aes_string(x = 'Date', y = factor))+
    geom_line()+
    xlab("")+
    ylab("")+
    ggtitle("Standardized residuals")+
    coord_cartesian(ylim = c(-5,5))+
    scale_y_continuous()
  
  out.newsimpact <- ggplot(newsimpdf, aes(x = x, y = y))+
    geom_line()+
    xlab(expression(epsilon[t - 1]))+
    ylab(expression(sigma[t]^2))+
    coord_cartesian(ylim = c(0,0.025))+
    ggtitle("News impact curve")
  
  
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
    ggtitle('Empirical density vs distributions')
  
  # Grid plots and print to pdf
  
  g <- arrangeGrob(out.ret, out.qq, out.acf, out.pacf,
                   out.aacf, out.apacf, out.newsimpact, out.empiricaldensity,
                   ncol = 2)
  ggsave(file= paste('output/garch_diagnostics/garch_diagnostics_', factor, '.jpeg', sep = ''), g, width = 8.3, height = 11.7, units = 'in', limitsize = F) #saves g
}