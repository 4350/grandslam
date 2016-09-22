#' Summary statistics plots
#' 
#' Creates 8 nice graphs of stats on returns and saves a jpeg with all
#' 1) Return series
#' 2) QQ-plot
#' 3) ACF
#' 4) PACF
#' 5) Absolute ACF
#' 6) Absolute PACF
#' 7) Volatility clustering (100 most extreme returns)
#' 8) Negative volatility clustering (100 most extreme neg returns)
#' 
#' @param df Data frame of returns
#' @param outlabel String indicating which data set used
#' @param factor Factor within the data frame to consider
#' 
#' @return g Grid 4x2 of the graphs
#' @export

statPlots <- function(df, outlabel, factor) {
  ret <- simplify2array(df[,factor])
  factorname <- factor
  absret <- abs(ret)
  
  out.acf <- autoplot(
    acf(ret, lag.max = 20, ylab = "", xlab = ""),
    main = "ACF log returns",
    ylab = "", 
    xlab = ""
  )
  
  out.pacf <- autoplot(
    pacf(ret, lag.max = 20, ylab = "", xlab = ""),
    main = "PACF log returns",
    ylab = "", 
    xlab = ""
  )
  
  out.aacf <- autoplot(
    acf(absret, lag.max = 20),
    main = "ACF absolute log returns",
    ylab = "", 
    xlab = ""
  )
  
  out.apacf <- autoplot(
    pacf(absret, lag.max = 20, ylab = "", xlab = ""),
    main = "PACF absolute log returns",
    ylab = "", 
    xlab = ""
  )
  
  out.qq <- ggplot(df, aes_string(sample = factor))+ 
    stat_qq(distribution = qnorm)+
    scale_y_continuous(labels = percent)+
    ggtitle("QQ plot vs normal distribution")
  
  out.ret <- ggplot(df, aes_string(x = 'Date', y = factor))+
    geom_line()+
    xlab("")+
    ylab("")+
    ggtitle("Log returns")+
    coord_cartesian(ylim = c(-.1,.1))+
    scale_y_continuous(labels = percent)
  
  # Volatility clustering
  extreme100limit <- tail(sort(absret), 100)[1]
  idDate <- absret >= extreme100limit
  df$extreme100 <- absret * idDate
  out.vol <- ggplot(df, aes(x = Date, y = extreme100))+
    geom_line()+
    xlab("")+
    ylab("")+
    coord_cartesian(ylim = c(0, 0.10))+
    scale_y_continuous(labels = percent)+
    ggtitle("100 most extreme returns")
  
  # Negative strings
  negextreme100limit <- head(sort(ret), 100)[100]
  idDate <- ret <= negextreme100limit
  df$negextreme100 <- ret * idDate
  out.negvol <- ggplot(df, aes(x = Date, y = negextreme100))+
    geom_line()+
    xlab("")+
    ylab("")+
    coord_cartesian(ylim = c(-.10, 0))+
    scale_y_continuous(labels = percent)+
    ggtitle("100 most extreme negative returns")
  
  # Grid plots and print to jpeg
  
  g <- arrangeGrob(out.ret, out.qq, out.acf, out.pacf, out.aacf, out.apacf, out.vol, out.negvol,
                   ncol = 2)
  ggsave(file= paste('output/MarginalStats/MarginalStats', factor, outlabel, 'jpeg', sep = '.'), g, width = 8.3, height = 11.7, units = 'in', limitsize = F) #saves g
  g
  
}