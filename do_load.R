#' Load daily factors return data and save weekly log returns
#' in two different data sets,
#' 1) full data set (1963-) weekly-full.RData
#' 2) estimation data set (1963-2010)`weekly-estim.RData

# Load Libraries ----
library(dplyr)

# Reset and Restrict Dataset ----
rm(list = ls())

to.weekly <- function(df) {
  log.return <- function(r) log(r / 100 + 1)
  
  df %>%
    mutate(
      # Convert Date into proper Date
      Date = as.Date(as.character(Date), format = "%Y%m%d"),
      # ...and get it as the Week
      Date = as.Date(cut(Date, "week")) + 4
    ) %>%
    
    # Get Weekly log total returns
    group_by(Date) %>% dplyr::summarise_each(funs(sum(log.return(.))))
}

df <-
  left_join(
    read.csv("data/source/ff_5factors-daily.csv") %>% to.weekly,
    read.csv("data/source/ff_momentum-daily.csv") %>% to.weekly,
    by = "Date"
  ) %>%
  select(-RF)

# Create two data sets, one for full period and one for estimation window (1963-2010)
df.estim <- df %>%
  filter(Date >= '1963-07-05' & Date <= '2010-12-31')

# Save for later loading
save(df, file = "data/derived/weekly-full.RData")
save(df.estim, file = "data/derived/weekly-estim.RData")
rm(to.weekly)
