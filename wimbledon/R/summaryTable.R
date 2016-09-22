#' Summary statistics table
#' 
#' Creates summary statistics table and writes to file
#' 
#' @param df Data frame of returns
#' 
#' @return out.table Summary table
#' @export

summaryTable <- function(df, outfilename) {
  out.table <- df %>%
    select(-Date) %>%
    basicStats() %>%
    .[c('nobs','Maximum','Minimum','Mean','Median','Stdev','Skewness','Kurtosis'),] %>%
    round(., digits = 4)
  
  write.table(out.table, file = paste('output/MarginalStats/summaryTable', outfilename , 'csv', sep = '.'), sep = ',')
  return(out.table)
}