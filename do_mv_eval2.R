rm(list = ls())

load('data/derived/distributions/oos_indep_10000.RData')
indep <- distribution
load('data/derived/distributions/oos_dynamic_std_10000.RData')
dstd <- distribution
rm(distribution)

plot(indep[ , 2, 900], indep[, 3, 900], ylim = c(-0.06, 0.06), xlim = c(-0.06, 0.06))
plot(dstd[ , 2, 900], dstd[, 3, 900], ylim = c(-0.06, 0.06), xlim = c(-0.06, 0.06))

ggplot(data.frame(x = indep[, 2, 900], y = indep[, 3, 900])) +
  geom_density2d(aes(x = x, y = y))+
  coord_cartesian(xlim = c(-.03,.03), ylim =c(-.03,.03))

ggplot(data.frame(x = dstd[, 2, 900], y = dstd[, 3, 900])) +
  geom_density2d(aes(x = x, y = y))+
  coord_cartesian(xlim = c(-.03,.03), ylim =c(-.03,.03))


# Correlations -----------------------------------------------------------

cors <- function(i, j) {
  sapply(seq(914), function(t) cor(dstd[, i, t], dstd[, j, t]))
}

plot(cors(1, 2)[394:914], type = 'l', ylim = c(-1, 1))
xplot(cors(2, 4)[576:914], type = 'l', ylim = c(-1, 1))
