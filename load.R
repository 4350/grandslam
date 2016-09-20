#' Load daily factors return data and save weekly log returns to
#'  `data/ff-weekly.RData`

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
  
  # Filter to Christoffersen & Langlois sample
  filter(Date >= '1963-07-05' & Date <= '2010-12-31')

# Save this for later loading
save(df, file = "data/ff-weekly.RData")
rm(to.weekly)

