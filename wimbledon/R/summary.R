#' Functionality related to summary statistics
#' 

#' Summary statistics table
#' 
#' Creates summary statistics table and writes to file
#' 
#' @param df Data frame of returns
#' 
#' @return out.table Summary table
#' @export

summary.table <- function(df, outfilename) {
  out.table <- df %>%
    select(-Date) %>%
    basicStats() %>%
    .[c('nobs','Maximum','Minimum','Mean','Median','Stdev','Skewness','Kurtosis'),] %>%
    round(., digits = 4)
  
  write.table(out.table, file = paste('output/MarginalStats/summaryTable', outfilename , 'csv', sep = '.'), sep = ',')
  return(out.table)
}

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

summary.plots <- function(df, outlabel, factor) {
  df <- as.data.frame(df)
  ret <- df[,factor]
  factorname <- factor
  absret <- abs(ret)
  
  out.acf <- autoplot(
    acf(ret, lag.max = 20, ylab = "", xlab = ""),
    main = "ACF log returns",
    ylab = "",
    xlab = ""
  ) +
    theme_Publication()
  
  out.pacf <- autoplot(
    pacf(ret, lag.max = 20, ylab = "", xlab = ""),
    main = "PACF log returns",
    ylab = "",
    ylim = c(-0.10,0.15),
    xlab = ""
  ) +
    theme_Publication()
  
  out.aacf <- autoplot(
    acf(absret, lag.max = 20),
    main = "ACF absolute log returns",
    ylab = "", 
    xlab = ""
  ) +
    theme_Publication()
  
  out.apacf <- autoplot(
    pacf(absret, lag.max = 20, ylab = "", xlab = ""),
    main = "PACF absolute log returns",
    ylab = "", 
    ylim = c(0,0.40),
    xlab = ""
  ) +
    theme_Publication()
  
  out.qq <- ggplot(df, aes_string(sample = factor))+ 
    stat_qq(distribution = qnorm)+
    scale_y_continuous(labels = percent)+
    ggtitle("QQ plot vs normal distribution")+
    theme_Publication()
  
  out.ret <- ggplot(df, aes_string(x = 'Date', y = factor))+
    geom_line()+
    xlab("")+
    ylab("")+
    ggtitle("Log returns")+
    coord_cartesian(ylim = c(-.1,.1))+
    scale_y_continuous(labels = percent)+
    theme_Publication()
  
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
    ggtitle("100 most extreme returns")+
    theme_Publication()
  
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
    ggtitle("100 most extreme negative returns")+
    theme_Publication()
  
  # Grid plots and print to jpeg
  
  g <- arrangeGrob(out.ret, out.qq, out.acf, out.pacf, out.aacf, out.apacf, out.vol, out.negvol,
                   ncol = 2)
  ggsave(file= paste('output/MarginalStats/MarginalStats', factor, outlabel, 'jpeg', sep = '.'), g, width = 14, height = 21, units = 'cm', limitsize = F) #saves g
  g
  
}


#' Cumulative return series plots
#' 
#' Creates 2 nice graphs of cumulative returns
#' 1) Cumulative returns to all strategies
#' 2) Cumulative standardized returns (10% annual vol) to all strategies 
#' Prints graphs to jpeg file
#' 
#' 
#' @param df Data frame of returns
#' @param outlabel String indicating which data set used
#' 
#' @return g Grid 2x1 of the graphs
#' @export
#' 
summary.cumretplots <- function(df, outlabel) {
  # Long format data series
  df <- gather(df, "factor", "ret", 2:ncol(df))
  
  # Mutate to create cumulative and standardized returns at 10% annual vol
  df <- df %>% 
    group_by(factor) %>%
    mutate(cumret = cumprod(1+ret)) %>%
    mutate(stdret = ret/sd(ret)*.10/sqrt(252)) %>%
    mutate(cumstdret = cumprod(1+stdret))
  
  # Plot cumulative return series
  out.cumret <- df %>%
    ggplot(aes(x=Date, y=cumret, group=factor)) +
    geom_line(aes(color=factor))+
    ggtitle('Cumulative returns to factor strategies')+
    xlab('')+
    ylab('Cumulative return')+
    theme(legend.position="bottom")+
    scale_x_date(date_breaks = "2 years", date_labels = "%y")+
    theme_Publication()
  
  # Plot cumulative standardized return series
  out.cumstdret <- df %>%
    ggplot(aes(x=Date, y=cumstdret, group=factor)) +
    geom_line(aes(color=factor))+
    ggtitle('Cumulative returns to factor strategies (standardized 10% annual vol)')+
    xlab('')+
    ylab('Cumulative return')+
    theme(legend.position="bottom")+
    scale_x_date(date_breaks = "2 years", date_labels = "%y")+
    theme_Publication()
  
  g <- arrangeGrob(out.cumret, out.cumstdret)
  ggsave(filename = paste('output/cumretPlot', outlabel, 'jpeg', sep = '.'), g, 'jpeg', dpi = 300, width = 14, height = 21, units = "cm", limitsize = F)
  g
  
}