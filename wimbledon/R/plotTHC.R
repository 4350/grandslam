#' Plot threshold correlation. To be used in list apply
#' 
#' @param value Data frame for one factor of threshold correlation data
#'
#' @return gridded plot of threshold correlation graphs
#' @export
plotTHC <- function(value) {
  ggplot(
    value,
    aes(x = qs, y = value, group = factor)
  ) + 
    geom_ribbon(aes(ymin = lb, ymax = ub),
                fill = "grey80") +
    geom_line(aes(color = factor)) +
    theme(legend.position="none") +
    ylab('') +
    coord_cartesian(xlim = c(0, 1), ylim = c(-0.5, 1)) +
    scale_x_continuous(labels = scales::percent) +
    facet_grid(. ~ order)  
}