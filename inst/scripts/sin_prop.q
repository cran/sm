n <- 50
x <- seq(0, 1, length = n)
m <- sin(2 * pi * x)
h <- 0.05
sigma <- 0.2
sm.sigma.old <- sm.sigma
sm.sigma <- function(x, y) { sigma }
model <- sm.regression(x, m, h = h, display = "none")
upper <- model$estimate + 2 * model$se
lower <- model$estimate - 2 * model$se
y <- rnorm(n, m, sigma)
plot(range(x), range(y, upper, lower), type = "n", 
	xlab="x", ylab="y")
polygon(c(x, rev(x)), c(upper, rev(lower)), border = F, col = "cyan")
sm.sigma <- sm.sigma.old
lines(x, m)
lines(x, model$estimate, lty = 3)
points(x, y)
