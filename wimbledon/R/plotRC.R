#' Plot rolling correlation. To be used in list apply
#' 
#' @param value Data frame for one factor of rolling correlation data
#'
#' @return gridded plot of rolling correlation graphs
#' @export
plotRC <- function(value) {
  ggplot(
    value,
    aes(x = Date, y = value, group = factor)
  ) + 
    geom_ribbon(aes(ymin = lb, ymax = ub),
                fill = "grey70") +
    geom_line(aes(color = factor)) +
    
    theme(legend.position="none") +
    ylab('') +
    coord_cartesian(ylim = c(-1, 1)) +
    #scale_x_date(date_breaks = "5 years", date_labels = "%y")+
    facet_grid(. ~ order)  
}