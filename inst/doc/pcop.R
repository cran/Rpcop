## -----------------------------------------------------------------------------
library(Rpcop)

set.seed(1)
n <- 120
t <- runif(n, -1, 1)
x <- cbind(t, t^2 + rnorm(n, sd = 0.08))

fit <- pcop(x, Ch = 1.5, Cd = 0.3, plot.true = FALSE)
names(fit)
head(fit$pcop.f1)
summary(fit)

## -----------------------------------------------------------------------------
str(fit$pcop.f2, max.level = 1)

## ----fig.width=5, fig.height=4------------------------------------------------
pcop(x, plot.true = TRUE, lwd = 2, col = 2)

