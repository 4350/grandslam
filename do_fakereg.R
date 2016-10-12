library(stargazer)

rm(list =ls())


# Fake garch table --------------------------------------------------------

mydata <- mtcars

#[1] mu     ar1     ma1   alpha1  beta1 eta11    skew shape omega 
#[1]"cyl"  "disp" "hp"   "drat"   "wt"  "qsec"  "vs"   "am" "gear"
colnames(mydata) <- c('y','mu', 'ar1','ma1','alpha1','beta1','eta11','skew','shape','omega')

Mkt.RF <- lm(y ~ 0 + mu + ar1 + ma1 + alpha1 + beta1 + eta11 + skew + shape + omega, data=mydata)
HML <- lm(y ~ 0 + mu + ar1 + ma1 + alpha1 + beta1 + eta11 + skew + shape + omega, data=mydata)
SMB <- lm(y ~ 0 + mu + ar1 + ma1 + alpha1 + beta1 + eta11 + skew + shape + omega, data=mydata)
Mom <- lm(y ~ 0 + mu + ar1 + ma1 + alpha1 + beta1 + eta11 + skew + shape + omega, data=mydata)
RMW <- lm(y ~ 0 + mu + ar1 + ma1 + alpha1 + beta1 + eta11 + skew + shape + omega, data=mydata)
CMA <- lm(y ~ 0 + mu + ar1 + ma1 + alpha1 + beta1 + eta11 + skew + shape + omega, data=mydata)

stargazer(Mkt.RF, HML, SMB,
      align = TRUE,
      omit.stat = c('adj.rsq','ser','f'),
      dep.var.caption = 'Factor series',
      object.names = TRUE,
      dep.var.labels.include = FALSE,
      notes.align = 'c',
      digits = 3,
      digit.separator = ','
)


# Fake copula table -------------------------------------------------------

colnames(mydata) <- c('y','df', 'gamma1','gamma2','gamma3','gamma4','gamma5','gamma6','alpha','beta')
cGauss <- lm(y ~ 0 + df + gamma1 + gamma2 + gamma3 + gamma4 + gamma5 + gamma6 + alpha + beta, data=mydata)
cGhst <- lm(y ~ 0 + df + gamma1 + gamma2 + gamma3 + gamma4 + gamma5 + gamma6 + alpha + beta, data=mydata)
cGhskt <- lm(y ~ 0 + df + gamma1 + gamma2 + gamma3 + gamma4 + gamma5 + gamma6 + alpha + beta, data=mydata)
dGauss <- lm(y ~ 0 + df + gamma1 + gamma2 + gamma3 + gamma4 + gamma5 + gamma6 + alpha + beta, data=mydata)
dGhst <- lm(y ~ 0 + df + gamma1 + gamma2 + gamma3 + gamma4 + gamma5 + gamma6 + alpha + beta, data=mydata)
dGhskt <- lm(y ~ 0 + df + gamma1 + gamma2 + gamma3 + gamma4 + gamma5 + gamma6 + alpha + beta, data=mydata)

stargazer(cGauss, cGhst, cGhskt,
          align = TRUE,
          omit.stat = c('adj.rsq','ser','f'),
          dep.var.caption = 'Copula specification',
          object.names = TRUE,
          dep.var.labels.include = FALSE,
          notes.align = 'c',
          digits = 3,
          digit.separator = ','
)


