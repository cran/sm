provide.data(tephra)
logit <- log(Al2O3/(100-Al2O3))
par(mfrow=c(1,2))
h.cv <- hcv(logit, display = "lines", ngrid = 32)
n  <- length(logit)
sd <- sqrt(var(logit))
h  <- seq(0.003, 0.054, length=32)
lines(h, nmise(sd, n, h) - 5.5, lty = 3)
sm.density(logit, h.cv)
sm.density(logit, lty = 3, add = T)
par(mfrow=c(1,1))
