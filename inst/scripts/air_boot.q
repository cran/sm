provide.data(aircraft)
y <- log(Span)[Period==3]
sm.density(y, xlab = "Log span")
for (i in 1:20) sm.density(sample(y, replace=T), col=6, add=T)
sm.density(y, xlab = "Log span", add=T)
