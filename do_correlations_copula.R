#' Copula correlations
#' 

library(devtools)
library(dplyr)
library(tidyr)
library(ggplot2)
load_all('wimbledon')

rm(list = ls())

MODEL_NAME <- 'dynamic_ghskt'
FACTORS <- c("Mkt.RF", "HML", "SMB", "Mom", "RMW", "CMA")

load(file.path('data/derived', sprintf('model_copula_%s.RData', MODEL_NAME)))
model.copula <- get(sprintf('model.copula.%s', gsub('_', '.', MODEL_NAME)))
colnames(model.copula$Correlation) <- FACTORS
rownames(model.copula$Correlation) <- FACTORS

load('data/derived/weekly-full.RData')

# Now we just get model.copula$Correlation[1, 2, ] for the correlation between
# factors 1 and 2. ez pz
correlation.df <- function(factor) {
  series <- t(model.copula$Correlation[factor,,])
  df <- data.frame(Date = df$Date, series)
  
  # Don't include correlation with ourselves
  df[[factor]] <- NA
  df
}

HML <- correlation.df('HML')
RMW <- correlation.df('RMW')
CMA <- correlation.df('CMA')

correlations <- bind_rows(HML = HML, RMW = RMW, CMA = CMA, .id = 'factor')
correlations$factor <- factor(correlations$factor, c("HML", "RMW", "CMA"))
gathered <- gather(correlations, factor2, value, 3:8, factor_key = TRUE)

# Get correlation version of Omega to represent unconditional expectation
unconditional <- data.frame(
  dc.Correlation(array(model.copula$Omega, dim = c(6, 6, 1)))[,, 1])

#unconditional[unconditional == 1] <- NA

colnames(unconditional) <- FACTORS
unconditional$factor2 <- factor(FACTORS)

unconditional <- unconditional[, c('HML', 'RMW', 'CMA', 'factor2')]

unconditional <- gather(unconditional, factor, value, 1:3, factor_key = TRUE)
unconditional$value <- round(unconditional$value, 2)

# Plotting Time ----------------------------------------------------------

g <- ggplot(gathered, aes(x = Date, y = value)) +
  geom_line(aes(color = 'Copula Correlation')) +
  geom_text(data = unconditional,
            aes(x = as.Date('2010-01-01'),
                y = -0.90,
                label = paste('r =', value)),
            family = 'Minion Pro', size = 3, parse = FALSE) +
  facet_grid(factor2 ~ factor) +
  theme_Publication() +
  scale_colour_Publication() +
  ylab('Correlation') +
  xlab('Year') + 
  scale_x_date(date_labels = "%y") +
  coord_cartesian(ylim = c(-1, 1), xlim = c(df$Date[1], df$Date[length(df$Date)]))

ggsave(sprintf('output/rollingCorrelations/copula_%s.png', MODEL_NAME),
       g, device = 'png', width = 14, height = 18, units = 'cm', dpi = 300, limitsize = F)