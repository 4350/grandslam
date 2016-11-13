#' Load daily factors return data and save weekly log returns
#' in two different data sets,
#' 1) full data set (1963-) weekly-full.RData
#' 2) estimation data set (1963-)`weekly-estim.RData

# Load Libraries ----
library(dplyr)

# Reset and Restrict Dataset ----
rm(list = ls())

to.weekly <- function(df) {
  gross.return <- function(r) (r / 100 + 1)
  
  df %>%
    mutate(
      # Convert Date into proper Date
      Date = as.Date(as.character(Date), format = "%Y%m%d"),
      # ...and get it as the Week
      Date = as.Date(cut(Date, "week")) + 4
    ) %>%
    
    # Get Weekly log total returns
    group_by(Date) %>% dplyr::summarise_each(funs(prod(gross.return(.)) - 1))
}

df <-
  left_join(
    read.csv("data/source/ff_5factors-daily.csv") %>% to.weekly,
    read.csv("data/source/ff_momentum-daily.csv") %>% to.weekly,
    by = "Date"
  ) %>%
  select(-RF)

# Explicitly order df as the order we use later
df <- select(
  df,
  Date,
  Mkt.RF,
  HML,
  SMB,
  Mom,
  RMW,
  CMA
)

# Create weekly data set
df.estim <- df %>%
  dplyr::filter(Date >= '1963-07-05')

# Get risk-free data set
df.RF <-
  left_join(
    read.csv("data/source/ff_5factors-daily.csv") %>% to.weekly,
    read.csv("data/source/ff_momentum-daily.csv") %>% to.weekly,
    by = "Date"
  ) %>%
  select(Date, RF) %>%
  filter(Date >= '1963-07-05')

# Save for later loading
save(df.estim, file = "data/derived/weekly-estim.RData")
save(df.RF, file = "data/derived/weekly-RF.RData")
rm(to.weekly)

# Load US recession data from NBER and convert to weekly data set ---------

usrec <- read.csv('data/source/fed_USRECD-daily.csv')
usrec <- usrec %>% mutate(DATE = as.Date(DATE), 
                          Date = as.Date(cut(DATE, 'week')) + 4,
                          recdummy = USRECD
                          ) %>%
  filter(DATE == Date) %>%
  select(Date, recdummy) %>%
  filter(Date >= '1963-07-05')

# Save for later loading
save(usrec, file = "data/derived/usrec-weekly.RData")
