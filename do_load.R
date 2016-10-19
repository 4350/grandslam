#' Load daily factors return data and save weekly log returns
#' in two different data sets,
#' 1) full data set (1963-) weekly-full.RData
#' 2) estimation data set (1963-)`weekly-estim.RData

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

# Create two weekly data sets, one for full period and one for estimation window (1963-)
df.estim <- df %>%
  dplyr::filter(Date >= '1963-07-05')

# Save for later loading
save(df, file = "data/derived/weekly-full.RData")
save(df.estim, file = "data/derived/weekly-estim.RData")
rm(to.weekly)

# Create two daily data sets, one for full period and one for estimation window (1963-)

df <-
  left_join(
    read.csv("data/source/ff_5factors-daily.csv"),
    read.csv("data/source/ff_momentum-daily.csv"),
    by = "Date"
  ) %>%
  select(-RF)

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

df <- df %>%
  mutate(
    # Convert Date into proper Date
    Date = as.Date(as.character(Date), format = "%Y%m%d")
  )

df.estim <- df %>%
  dplyr::filter(Date >= '1963-07-05')


# Save for later loading
save(df, file = "data/derived/daily-full.RData")
save(df.estim, file = "data/derived/daily-estim.RData")


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
