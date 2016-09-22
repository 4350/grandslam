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
cumretPlots <- function(df, outlabel) {
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
    geom_line(aes(linetype=factor, color=factor))+
    geom_point(aes(color=factor))+
    ggtitle('Cumulative returns to factor strategies')+
    xlab('')+
    ylab('Cumulative return')+
    theme(legend.position="bottom")+
    scale_x_date(date_breaks = "2 years", date_labels = "%y")
  
  # Plot cumulative standardized return series
  out.cumstdret <- df %>%
    ggplot(aes(x=Date, y=cumstdret, group=factor)) +
    geom_line(aes(linetype=factor, color=factor))+
    geom_point(aes(color=factor))+
    ggtitle('Cumulative returns to factor strategies (standardized 10% annual vol)')+
    xlab('')+
    ylab('Cumulative return')+
    theme(legend.position="bottom")+
    scale_x_date(date_breaks = "2 years", date_labels = "%y")
  
  g <- arrangeGrob(out.cumret, out.cumstdret)
  ggsave(filename = paste('output/cumretPlot', outlabel, 'jpeg', sep = '.'), g, 'jpeg', dpi = 300, width = 8.3, height = 11.7, units = "in", limitsize = F)
  g
  
}