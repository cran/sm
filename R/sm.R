"addplot" <-
function (d, f, theta, phi)
{
    a <- (f * pi)/180
    b <- (d * pi)/180
    radtheta <- (theta * pi)/180
    radphi <- (phi * pi)/180
    xyzcheck <<- ((cos(a) * cos(b) * sin(radtheta) * cos(radphi)) +
        (sin(b) * sin(radphi)) - (sin(a) * cos(b) * cos(radphi) *
        cos(radtheta)))
    llong <<- a[xyzcheck >= 0]
    llat <<- b[xyzcheck >= 0]
    if (length(llat) == 0) {
        break
    }
    X <<- (cos(llong) * cos(llat) * cos(radtheta)) + (sin(llong) *
        cos(llat) * sin(radtheta))
    Y <<- (sin(llat) * cos(radphi)) + ((cos(llat) * sin(radphi)) *
        ((sin(llong) * cos(radtheta)) - (cos(llong) * sin(radtheta))))
    par(pty = "s")
    points(X, Y)
}
"ask" <-
function (message = "Type in datum")
eval(parse(prompt = paste(message, ": ", sep = "")))
"attach.frame" <-
function (data, name, ...)
{
    if (missing(name))
        name <- deparse(substitute(data))
    if (is.data.frame(data)) {
        if (!is.na(pos <- match(name, search()))) {
            cat(paste(name, "already attached, re-attached in 2nd position\n"))
            detach(pos = pos)
        }
        cat(paste("attaching", name, "\n", sep = " "))
        attach(what = data, pos = 2, name = name, ...)
    }
    else {
        cat(name)
        cat(" is not a data.frame\n")
    }
    invisible()
}
"binning" <-
function (x, y, breaks, nbins)
{
    binning.1d <- function(x, y, breaks, nbins) {
        f <- cut(x, breaks = breaks)
        if (any(is.na(f)))
            stop("breaks do not span the range of x")
        freq <- tabulate(f, length(levels(f)))
        midpoints <- (breaks[-1] + breaks[-(nbins + 1)])/2
        id <- (freq > 0)
        x <- midpoints[id]
        x.freq <- as.vector(freq[id])
        result <- list(x = x, x.freq = x.freq, table.freq = freq,
            breaks = breaks)
        if (!all(is.na(y))) {
            result$means <- as.vector(tapply(y, f, mean))[id]
            result$sums <- as.vector(tapply(y, f, sum))[id]
            result$devs <- as.vector(tapply(y, f, function(x) sum((x -
                mean(x))^2)))[id]
        }
        result
    }
    binning.2d <- function(x, y, breaks, nbins) {
        f1 <- cut(x[, 1], breaks = breaks[, 1])
        f2 <- cut(x[, 2], breaks = breaks[, 2])
        freq <- t(table(f1, f2))
        dimnames(freq) <- NULL
        midpoints <- (breaks[-1, ] + breaks[-(nbins + 1), ])/2
        z1 <- midpoints[, 1]
        z2 <- midpoints[, 2]
        X <- cbind(rep(z1, length(z2)), rep(z2, rep(length(z1), length(z2))))
        X.f <- as.vector(t(freq))
        id <- (X.f > 0)
        X <- X[id, ]
        dimnames(X) <- list(NULL, dimnames(x)[[2]])
        X.f <- X.f[id]
        result <- list(x = X, x.freq = X.f, midpoints = midpoints,
            breaks = breaks, table.freq = freq)
        if (!all(is.na(y))) {
            result$means <- as.numeric(tapply(y, list(f1, f2),
                mean))[id]
            result$devs <- as.numeric(tapply(y, list(f1, f2),
                function(x) sum((x - mean(x))^2)))[id]
        }
        result
    }
    if (length(dim(x)) > 0) {
        if (!isMatrix(x))
            stop("wrong parameter x for binning")
        if (dim(x)[2] != 2)
            stop("wrong parameter x for binning")
        if (missing(y))
            y <- rep(NA, nrow(x))
        if (missing(nbins))
            nbins <- round(log(nrow(x))/log(2) + 1)
        if (missing(breaks)) {
            breaks <- cbind(seq(min(x[, 1]), max(x[, 1]), length = nbins + 1),
                            seq(min(x[, 2]), max(x[, 2]), length = nbins + 1))
            breaks[1, ] <- breaks[1, ] - rep(10^(-5), 2)
        }
        else nbins <- length(breaks)/2 - 1
        if (max(abs(breaks)) == Inf | is.na(max(abs(breaks))))
            stop("Illegal breaks")
        result <- binning.2d(x, y, breaks = breaks, nbins = nbins)
    }
    else {
        x <- as.vector(x)
        if (missing(y))
            y <- rep(NA, length(x))
        if (missing(nbins))
            nbins <- round(log(length(x))/log(2) + 1)
        if (missing(breaks)) {
            breaks <- seq(min(x), max(x), length = nbins + 1)
            breaks[1] <- breaks[1] - 10^(-5)
        }
        else nbins <- length(breaks) - 1
        if (max(abs(breaks)) == Inf | is.na(max(abs(breaks))))
            stop("Illegal breaks")
        result <- binning.1d(x, y, breaks = breaks, nbins = nbins)
    }
    result
}
"box1." <-
function (theta = pi/6, phi = pi/6, sc, col = par("col"), axes.lim = sc)
{
    plot(c(-sqrt(3), sqrt(3)), c(-sqrt(3), sqrt(3)), type = "n",
        xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "")
    lines(rbind(coord2.(axes.lim[1, 1], axes.lim[2, 1], axes.lim[3,
        1], theta, phi, sc), coord2.(axes.lim[1, 1], axes.lim[2,
        1], axes.lim[3, 2], theta, phi, sc), coord2.(axes.lim[1,
        2], axes.lim[2, 1], axes.lim[3, 2], theta, phi, sc)),
        lty = 3)
    lines(rbind(coord2.(axes.lim[1, 1], axes.lim[2, 1], axes.lim[3,
        2], theta, phi, sc), coord2.(axes.lim[1, 1], axes.lim[2,
        2], axes.lim[3, 2], theta, phi, sc)), lty = 3)
}
"box2." <-
function (theta = pi/6, phi = pi/6, sc, labels = c("", "", ""),
    col = par("col"), cex = 9/10, axes.lim = sc)
{
    lines(rbind(coord2.(axes.lim[1, 2], axes.lim[2, 1], axes.lim[3,
        1], theta, phi, sc), coord2.(axes.lim[1, 2], axes.lim[2,
        2], axes.lim[3, 1], theta, phi, sc), coord2.(axes.lim[1,
        1], axes.lim[2, 2], axes.lim[3, 1], theta, phi, sc),
        coord2.(axes.lim[1, 1], axes.lim[2, 1], axes.lim[3, 1],
            theta, phi, sc), coord2.(axes.lim[1, 2], axes.lim[2,
            1], axes.lim[3, 1], theta, phi, sc), coord2.(axes.lim[1,
            2], axes.lim[2, 1], axes.lim[3, 2], theta, phi, sc),
        coord2.(axes.lim[1, 2], axes.lim[2, 2], axes.lim[3, 2],
            theta, phi, sc), coord2.(axes.lim[1, 1], axes.lim[2,
            2], axes.lim[3, 2], theta, phi, sc), coord2.(axes.lim[1,
            1], axes.lim[2, 2], axes.lim[3, 1], theta, phi, sc)),
        lty = 3)
    lines(rbind(coord2.(axes.lim[1, 2], axes.lim[2, 2], axes.lim[3,
        1], theta, phi, sc), coord2.(axes.lim[1, 2], axes.lim[2,
        2], axes.lim[3, 2], theta, phi, sc)), lty = 3)
    cex <- 0.9
    text(rbind(coord2.((axes.lim[1, 1] + axes.lim[1, 2])/2, axes.lim[2,
        1], axes.lim[3, 1], theta, phi, sc)), labels[1], adj = 1,
        cex = cex, col = col)
    text(rbind(coord2.(axes.lim[1, 1], (axes.lim[2, 1] + axes.lim[2,
        2])/2, axes.lim[3, 1], theta, phi, sc)), labels[2], adj = 1,
        cex = cex, col = col)
    text(rbind(coord2.(axes.lim[1, 2], axes.lim[2, 1], (axes.lim[3,
        1] + axes.lim[3, 2])/2, theta, phi, sc)), labels[3],
        adj = 0, cex = cex, col = col)
}
"britmap" <-
function ()
{
    provide.data(britpts)
    jump <- c(0, sqrt(diff(britlat)^2 + diff(britlong)^2))
    flag <- rep(1, nrow(britpts))
    flag[jump >= 0.6] <- NA
    lines(britpts * flag)
}
"change" <-
function (th, ph)
{
    cat("Theta =", th, "\n")
    cat("Phi =", ph, "\n")
    scan(n = 1)
    cat("Change theta to ? \n")
    theta <<- scan(n = 1)
    if (theta >= 360) theta <<- theta - 360
    cat("\n", "Change phi to ? \n")
    phi <<- scan(n = 1)
    if (phi > 90) phi <<- 90
    if (phi < -90) phi <<- -90
    cat("Theta =", theta, "\n")
    cat("Phi =", phi, "\n")
    list(theta = theta, phi = phi)
}
"circle" <-
function (r)
{
    angle <- seq(0, 7, by = 0.1)
    x <- r * cos(angle)
    y <- r * sin(angle)
    par(lty = 1)
    lines(x, y)
}
"coord2." <-
function (x, y, z, theta, phi, sc)
{
    x <- -1 + ((x - sc[1, 1]) * 2)/(sc[1, 2] - sc[1, 1])
    y <- -1 + ((y - sc[2, 1]) * 2)/(sc[2, 2] - sc[2, 1])
    z <- -((-1 + ((z - sc[3, 1]) * 2)/(sc[3, 2] - sc[3,
        1])))
    co <- cbind(x * cos(theta) - z * sin(theta), y * cos(phi) -
        (x * sin(theta) + z * cos(theta)) * sin(phi))
    co
}
"cv" <-
function (x, h, ...)
{
    opt <- sm.options(list(...))
    if (!isMatrix(x)) {
        n <- length(x)
        replace.na(opt, h.weights, rep(1, n))
        hcvff <- sum(dnorm(0, mean = 0, sd = sqrt(2) * h * opt$h.weights))/(n *
            (n - 1))
        W <- matrix(rep(x, rep(n, n)), ncol = n, byrow = TRUE)
        W <- W - matrix(rep(x, n), ncol = n, byrow = TRUE)
        W1 <- matrix(rep(opt$h.weights^2, n), ncol = n, byrow = TRUE)
        W2 <- exp(-0.5 * (W/(h * sqrt(W1 + t(W1))))^2)/(sqrt(2 * pi) *
                                                        h * sqrt(W1 + t(W1)))
        hcvff <- hcvff + (sum(W2) - sum(diag(W2))) * (n - 2)/(n * (n - 1)^2)
        W2 <- exp(-0.5 * (W/(h * sqrt(W1)))^2)/(sqrt(2 * pi) * h * sqrt(W1))
        hcvff <- hcvff - (sum(W2) - sum(diag(W2))) * 2/(n * (n - 1))
    }
    if (isMatrix(x)) {
        x1 <- x[, 1]
        x2 <- x[, 2]
        h1 <- h * sqrt(var(x1))
        h2 <- h * sqrt(var(x2))
        n <- length(x1)
        replace.na(opt, h.weights, rep(1, n))
        hcvff <- sum(dnorm(0, mean = 0, sd = sqrt(2) * h1 * opt$h.weights) *
            dnorm(0, mean = 0, sd = sqrt(2) * h2 * opt$h.weights))/(n *
            (n - 1))
        W <- matrix(rep(x1, rep(n, n)), ncol = n, byrow = TRUE)
        W <- W - matrix(rep(x1, n), ncol = n, byrow = TRUE)
        W1 <- matrix(rep(opt$h.weights^2, n), ncol = n, byrow = TRUE)
        W2 <- exp(-0.5 * (W/(h1 * sqrt(W1 + t(W1))))^2)/(sqrt(2 *
            pi) * h1 * sqrt(W1 + t(W1)))
        W <- matrix(rep(x2, rep(n, n)), ncol = n, byrow = TRUE)
        W <- W - matrix(rep(x2, n), ncol = n, byrow = TRUE)
        W2 <- W2 * exp(-0.5 * (W/(h2 * sqrt(W1 + t(W1))))^2)/(sqrt(2 *
            pi) * h2 * sqrt(W1 + t(W1)))
        hcvff <- hcvff + (sum(W2) - sum(diag(W2))) * (n - 2)/(n *
            (n - 1)^2)
        W2 <- exp(-0.5 * (W/(h2 * sqrt(W1)))^2)/(sqrt(2 * pi) *
            h2 * sqrt(W1))
        W <- matrix(rep(x1, rep(n, n)), ncol = n, byrow = TRUE)
        W <- W - matrix(rep(x1, n), ncol = n, byrow = TRUE)
        W2 <- W2 * exp(-0.5 * (W/(h1 * sqrt(W1)))^2)/(sqrt(2 *
            pi) * h1 * sqrt(W1))
        hcvff <- hcvff - (sum(W2) - sum(diag(W2))) * 2/(n * (n - 1))
    }
    hcvff
}
"hcv" <-
function (x, y = NA, hstart = NA, hend = NA, ...)
{
    opt <- sm.options(list(...))
    replace.na(opt, ngrid, 8)
    replace.na(opt, display, "none")
    if (length(dim(x)) > 0) {
        ndim <- 2
        n <- length(x[, 1])
    }
    else {
        ndim <- 1
        n <- length(x)
    }
    replace.na(opt, h.weights, rep(1, n))
    ngrid <- opt$ngrid
    display <- opt$display
    h.weights <- opt$h.weights
    if (is.na(hstart)) {
        if (ndim == 1) hstart <- hnorm(x)/10
        else if (any(is.na(y)))
            hstart <- hnorm(x[, 1]/sqrt(var(x[, 1])))/10
        else hstart <- hnorm(x[, 1]/sqrt(var(x[, 1])))/4
    }
    if (is.na(hend)) {
        if (ndim == 1)
            hend <- hnorm(x) * 2
        else hend <- hnorm(x[, 1]/sqrt(var(x[, 1]))) * 2
    }
    cvgrid <- vector("numeric", length = ngrid)
    hgrid <- log(hstart) + (log(hend) - log(hstart)) * (0:(ngrid - 1))/
        (ngrid - 1)
    if (any(is.na(y))) {
        for (i in 1:ngrid) cvgrid[i] <- cv(x, exp(hgrid[i]),
            h.weights = h.weights)
    }
    else {
        if (ndim == 1)
            for (i in 1:ngrid) {
                cvgrid[i] <- sum((y - sm.weight(x, x, h = exp(hgrid[i]),
                  cross = TRUE, options = list(h.weights = h.weights)) %*% y)^2)
            }
        if (ndim == 2)
            for (i in 1:ngrid) {
                cvgrid[i] <-
                    sum((y - sm.weight2(x, x, exp(hgrid[i] * c(sqrt(var(x[, 1])), sqrt(var(x[, 2])))), h.weights = h.weights, cross = TRUE) %*% y)^2)
            }
    }
    if (any(is.na(cvgrid))) {
        cat("\n")
        cat("hcv: some computations failed.", "\n")
        cat("Try readjusting hstart and hend.", "\n")
        cat("hstart: ", hstart, "\n")
        cat("hend  : ", hend, "\n")
        cat("\n")
        print(cbind(h = exp(hgrid), cv = cvgrid))
        stop()
    }
    ind <- (1:ngrid)[cvgrid == min(cvgrid)]
    if (!(display == "none")) {
        if (!opt$add) {
            if (display == "log")
                plot(hgrid, cvgrid, type = "l", xlab = "Log h", ylab = "CV")
            else plot(exp(hgrid), cvgrid, type = "l", xlab = "h", ylab = "CV")
        }
        else {
            if (display == "log")
                lines(hgrid, cvgrid)
            else lines(exp(hgrid), cvgrid)
        }
    }
    if (ind == 1 | ind == ngrid) {
        cat("\n")
        cat("hcv: boundary of search area reached.", "\n")
        cat("Try readjusting hstart and hend.", "\n")
        cat("hstart: ", hstart, "\n")
        cat("hend  : ", hend, "\n")
        cat("\n")
        print(cbind(h = exp(hgrid), cv = cvgrid))
        stop()
    }
    v0 <- cvgrid[ind - 1]
    v1 <- cvgrid[ind]
    v2 <- cvgrid[ind + 1]
    l0 <- hgrid[ind - 1]
    l1 <- hgrid[ind]
    l2 <- hgrid[ind + 1]
    aa <- (v1 - v0 - (l1 - l0) * (v1 - v2)/(l1 - l2))/(l1^2 -
        l0^2 - (l1^2 - l2^2) * (l1 - l0)/(l1 - l2))
    bb <- (v1 - v2 - aa * (l1^2 - l2^2))/(l1 - l2)
    cc <- v0 - aa * l0^2 - bb * l0
    h <- exp(-bb/(2 * aa))
    if (ndim == 1) result <- h
    else result <- c(h * sqrt(var(x[, 1])), h * sqrt(var(x[, 2])))
    result
}
"hidplot" <-
function (invis, theta, phi)
{
    invislong <- invis$invislong
    invislat <- invis$invislat
    par(pch = "O")
    a <- (invislong * pi)/180
    b <- (invislat * pi)/180
    radtheta <- (theta * pi)/180
    radphi <- (phi * pi)/180
    if (length(invislat) == 0) {
        points(0, 0, type = "n")
        break
    }
    X <<- (cos(invislong) * cos(invislat) * cos(radtheta)) +
        (sin(invislong) * cos(invislat) * sin(radtheta))
    Y <<- (sin(invislat) * cos(radphi)) + ((cos(invislat) * sin(radphi)) *
        ((sin(invislong) * cos(radtheta)) - (cos(invislong) *
            sin(radtheta))))
    points(X, Y)
}
"hnorm" <-
function (x, weights)
{
    if (isMatrix(x)) {
        if (missing(weights))
            weights <- rep(1, nrow(x))
        ndim <- ncol(x)
        n <- sum(weights)
        sd <- sqrt(apply(x, 2, wvar, w = weights))
        if (ndim == 2)
            hh <- sd * (1/n)^(1/6)
        if (ndim == 3)
            hh <- sd * (4/(5 * n))^(1/7)
        hh
    }
    else {
        if (missing(weights))
            weights <- rep(1, length(x))
        sd <- sqrt(wvar(x, weights))
        sd * (4/(3 * sum(weights)))^(1/5)
    }
}
"hsj" <-
function (x)
{
    h0 <- hnorm(x)
    v0 <- sj(x, h0)
    if (v0 > 0) hstep <- 1.1
    else hstep <- 0.9
    h1 <- h0 * hstep
    v1 <- sj(x, h1)
    while (v1 * v0 > 0) {
        h0 <- h1
        v0 <- v1
        h1 <- h0 * hstep
        v1 <- sj(x, h1)
    }
    h0 + (h1 - h0) * abs(v0)/(abs(v0) + abs(v1))
}
"incphi" <-
function (ph, inc)
{
    phi <<- ph + inc
    if (phi > 90) phi <<- 90
    if (phi < -90) phi <<- -90
    cat("Phi =", phi, "\n")
    phi
}
"inctheta" <-
function (th, inc)
{
    theta <<- th + inc
    if (theta >= 360) theta <<- theta - 360
    theta
}
"latlines" <-
function (beta, theta, phi)
{
    if (beta < (phi - 90) | beta > (phi + 90)) return()
    par(pch = ".")
    radtheta <- (theta * pi)/180
    radbeta <- (beta * pi)/180
    radphi <- (phi * pi)/180
    alpha <- seq(0, (2 * pi), by = 0.05)
    xyzcheck <<- ((cos(alpha) * cos(radbeta) * sin(radtheta) *
        cos(radphi)) + (sin(radbeta) * sin(radphi)) - (sin(alpha) *
        cos(radbeta) * cos(radphi) * cos(radtheta)))
    alphaplot <- alpha[xyzcheck >= 0]
    X <- (cos(alphaplot) * cos(radbeta) * cos(radtheta)) + (sin(alphaplot) *
        cos(radbeta) * sin(radtheta))
    Y <- (sin(radbeta) * cos(radphi)) + (((sin(alphaplot) * cos(radtheta)) -
        (cos(alphaplot) * sin(radtheta))) * cos(radbeta) * sin(radphi))
    points(X, Y)
}
"latlines.e" <-
function (beta, theta, phi)
{
    if (beta < (phi - 90) | beta > (phi + 90)) return()
    par(lty = 2)
    par(pch = ".")
    radtheta <- (theta * pi)/180
    radbeta <- (beta * pi)/180
    radphi <- (phi * pi)/180
    alpha <- seq(0, (2 * pi), by = 0.005)
    xyzcheck <<- ((cos(alpha) * cos(radbeta) * sin(radtheta) *
        cos(radphi)) + (sin(radbeta) * sin(radphi)) - (sin(alpha) *
        cos(radbeta) * cos(radphi) * cos(radtheta)))
    alphaplot <- alpha[xyzcheck >= 0]
    X <- (cos(alphaplot) * cos(radbeta) * cos(radtheta)) + (sin(alphaplot) *
        cos(radbeta) * sin(radtheta))
    Y <- (sin(radbeta) * cos(radphi)) + (((sin(alphaplot) * cos(radtheta)) -
        (cos(alphaplot) * sin(radtheta))) * cos(radbeta) * sin(radphi))
    points(X, Y)
}
"longlines" <-
function (alpha, theta, phi)
{
    par(pch = ".")
    radtheta <- (theta * pi)/180
    radalpha <- (alpha * pi)/180
    radphi <- (phi * pi)/180
    beta <- seq(0, (2 * pi), by = 0.05)
    xyzcheck <<- ((cos(radalpha) * cos(beta) * sin(radtheta) *
        cos(radphi)) + (sin(beta) * sin(radphi)) - (sin(radalpha) *
        cos(beta) * cos(radphi) * cos(radtheta)))
    betaplot <- beta[xyzcheck >= 0]
    X <- (cos(radalpha) * cos(betaplot) * cos(radtheta)) + (sin(radalpha) *
        cos(betaplot) * sin(radtheta))
    Y <- (sin(betaplot) * cos(radphi)) + (((sin(radalpha) * cos(radtheta)) -
        (cos(radalpha) * sin(radtheta))) * cos(betaplot) * sin(radphi))
    points(X, Y)
}
"longlines.e" <-
function (alpha, theta, phi)
{
    par(lty = 2)
    par(pch = ".")
    radtheta <- (theta * pi)/180
    radalpha <- (alpha * pi)/180
    radphi <- (phi * pi)/180
    beta <- seq(0, (2 * pi), by = 0.005)
    xyzcheck <<- ((cos(radalpha) * cos(beta) * sin(radtheta) *
        cos(radphi)) + (sin(beta) * sin(radphi)) - (sin(radalpha) *
        cos(beta) * cos(radphi) * cos(radtheta)))
    betaplot <- beta[xyzcheck >= 0]
    X <- (cos(radalpha) * cos(betaplot) * cos(radtheta)) + (sin(radalpha) *
        cos(betaplot) * sin(radtheta))
    Y <- (sin(betaplot) * cos(radphi)) + (((sin(radalpha) * cos(radtheta)) -
        (cos(radalpha) * sin(radtheta))) * cos(betaplot) * sin(radphi))
    points(X, Y)
}
"nise" <-
function (y, ...)
{
    n <- length(y)
    opt <- sm.options(list(...))
    replace.na(opt, nbins, round((n > 500) * 8 * log(n)))
    replace.na(opt, hmult, 1)
    if (opt$nbins > 0) {
        bins <- binning(y, nbins = nbins)
        y <- bins$x
        weights <- bins$x.freq
    }
    else weights <- rep(1, n)
    y <- (y - wmean(y, weights))/sqrt(wvar(y, weights))
    h <- hnorm(y) * opt$hmult
    result <- dnorm(0, sd = sqrt(2 + 2 * h^2))
    result <- result - 2 * sm.density(y, h = sqrt(1 + 2 * h^2),
        eval.points = 0, display = "none", weights = weights,
        nbins = 0)$estimate
    result <- result + wmean(sm.density(y, h = sqrt(2) * h, eval.points = y,
        display = "none", weights = weights, nbins = 0)$estimate,
        weights)
    result
}
"nmise" <-
function (sd, n, h)
{
    dnorm(0, sd = sqrt(2) * h)/n + (1 - 1/n) * dnorm(0, sd = sqrt(2 *
        (sd^2 + h^2))) - 2 * dnorm(0, sd = sqrt(2 * sd^2 + h^2)) +
        dnorm(0, sd = sqrt(2) * sd)
}
"nnbr" <-
function (x, k)
{
    if (isMatrix(x)) {
        ndim <- 2
        n <- nrow(x)
    }
    else {
        ndim <- 1
        n <- length(x)
    }
    knn <- vector("numeric", n)
    if (ndim == 1) {
        for (i in 1:length(x)) knn[i] <- sort(abs(x - x[i]))[k +
            1]}
    if (ndim == 2) {
        for (i in 1:length(x[, 1])) knn[i] <- sort(sqrt(((x[,
            1] - x[i, 1])^2)/var(x[, 1]) + ((x[, 2] - x[i, 2])^2)/var(x[,
            2])))[k + 1]
    }
    knn
}
"normdens.band" <-
function (x, h, weights = rep(1, length(x)), options = list())
{
    opt <- sm.options(options)
    xlim <- opt$xlim
    yht <- opt$yht
    ngrid <- opt$ngrid
    x.points <- seq(xlim[1], xlim[2], length = ngrid)
    xbar <- wmean(x, weights)
    sx <- sqrt(wvar(x, weights))
    hm <- h * opt$hmult
    dmean <- dnorm(x.points, xbar, sqrt(sx^2 + hm^2))
    dvar <- (dnorm(0, 0, sqrt(2 * hm^2)) * dnorm(x.points, xbar,
        sqrt(sx^2 + 0.5 * hm^2)) - (dmean)^2)/sum(weights)
    upper <- pmin(dmean + 2 * sqrt(dvar), par()$usr[4])
    lower <- pmax(0, dmean - 2 * sqrt(dvar))
    polygon(c(par()$usr[1:2], par()$usr[2:1]), rep(c(par()$usr[3],
        par()$usr[4] * 0.999), c(2, 2)), col = 0, border= 0)
    polygon(c(x.points, rev(x.points)), c(upper, rev(lower)),
            col = "cyan", border = 0)
}
"np.contour.plot.3d." <-
function (coord, data = matrix(), shadow = TRUE, gridsize = 20,
    numpts = 3, xmin = NA, xmax = NA, ymin = NA, ymax = NA, zmin = NA,
    zmax = NA, xlab = NA, ylab = NA, zlab = NA, theta = pi/4,
    phi = pi/4, colour = FALSE, title.colour = 3, label.colour = 3,
    axes.colour = 6, plot.colour = "blue", shadow.colour = "cyan",
    cex = NA)
{
    col <- par("col")
    on.exit(par(col = col))
    if (!colour) {
        title.colour <- col
        label.colour <- col
        axes.colour <- col
    }
    opar <- par(pty = "s")
    on.exit(par(opar))
    if (is.na(xlab))
        if (!is.null(attributes(data)$dimnames)) {
            xlab <- attributes(data)$dimnames[[2]][1]
            ylab <- attributes(data)$dimnames[[2]][2]
            zlab <- attributes(data)$dimnames[[2]][3]
        }
        else {
            xlab <- ""
            ylab <- ""
            zlab <- ""
        }
    if (is.na(xmin)) {
        xmin <- min(data[, 1])
        xmax <- max(data[, 1])
        ymin <- min(data[, 2])
        ymax <- max(data[, 2])
        zmin <- min(data[, 3])
        zmax <- max(data[, 3])
    }
    sc <- matrix(c(xmin, xmax, ymin, ymax, zmin, zmax), ncol = 2,
        byrow = TRUE)
    dst <- sqrt(coord[, 1]^2 + coord[, 2]^2 + (zmax - coord[,
        3])^2)
    coord[, 1] <- xmin + (coord[, 1] - 1)/(gridsize - 1) * (xmax - xmin)
    coord[, 2] <- ymin + (coord[, 2] - 1)/(gridsize - 1) * (ymax - ymin)
    coord[, 3] <- zmin + (coord[, 3] - 1)/(gridsize - 1) * (zmax - zmin)
    poly <- matrix(t(coord), ncol = 3 * numpts, byrow = TRUE)
    numpoly <- nrow(poly)
    dst <- matrix(dst, ncol = numpts, byrow = TRUE)
    min.dst <- apply(dst, 1, min)
    xyz <- poly[order(min.dst), ]
    par(col = axes.colour)
    box1.(theta, phi, sc, title.colour)
    par(col = shadow.colour)
    if (shadow) {
        for (i in (1:numpoly)) {
            polyi <- matrix(xyz[i, ], ncol = 3, byrow = TRUE)
            polygon(coord2.(polyi[, 1], ymin, polyi[, 3], theta,
                phi, sc), col = shadow.colour, border = shadow.colour)
        }
    }
    par(col = plot.colour)
    for (i in (1:numpoly)) {
        polyi <- matrix(xyz[i, ], ncol = 3, byrow = TRUE)
        polygon(coord2.(polyi[, 1], polyi[, 2], polyi[, 3], theta, phi, sc),
                col = 0)
        polygon(coord2.(polyi[, 1], polyi[, 2], polyi[, 3], theta, phi, sc))
    }
    par(col = axes.colour)
    box2.(theta, phi, sc, c(xlab, ylab, zlab), label.colour, cex)
}
"p.quad.moment" <-
function (A, Sigma, tobs, ndevs)
{
    B <- A %*% Sigma
    k1 <- sum(diag(B)) - tobs * ndevs
    C <- B %*% B
    k2 <- 2 * sum(diag(C)) + 2 * tobs^2 * ndevs
    k3 <- 8 * sum(diag(C %*% B)) - 8 * tobs^3 * ndevs
    aa <- abs(k3/(4 * k2))
    bb <- (8 * k2^3)/k3^2
    cc <- k1 - aa * bb
    1 - pchisq(-cc/aa, bb)
}
"pause" <-
function ()
{
    cat("Pause. Press <Enter> to continue...")
    readline()
    invisible()
}
"smplot.density" <-
function (x, h, weights = rep(1, length(x)), rawdata = list(x = x),
    options = list())
{
    opt <- sm.options(options)
    if (opt$positive)
        est <- sm.density.positive.1d(x, h, weights = weights,
            options = opt)
    else {
        est <- sm.density.eval.1d(x, h, weights = weights, options = opt)
        if (opt$band)
            normdens.band(x, h, weights, options = opt)
        else if (!opt$add)
            polygon(c(par()$usr[1:2], par()$usr[2:1]), rep(c(par()$usr[3],
                par()$usr[4] * 0.999), c(2, 2)), col = 0, border = 0)
    }
    box()
    lines(est$eval.points, est$estimate, lty = opt$lty, col = opt$col)
    if (opt$rugplot && !opt$add)
        rug(jitter(rawdata$x, amount = 0), 0.015)
    if (opt$display %in% "se" & all(opt$h.weights == rep(1, length(x)))) {
        se <- sqrt(dnorm(0, sd = sqrt(2))/(4 * h * sum(weights)))
        upper <- sqrt(est$estimate) + 2 * se
        lower <- pmax(sqrt(est$estimate) - 2 * se, 0)
        upper <- upper^2
        lower <- lower^2
        lines(est$eval.points, upper, lty = 3, col = opt$col)
        lines(est$eval.points, lower, lty = 3, col = opt$col)
    }
    invisible(est)
}
"plot.regression" <-
function (x, y, design.mat, h, r, model, weights, rawdata = list(),
    options = list(), ...)
{
    opt <- sm.options(options)
    rnew <- sm.regression.eval.1d(x, y, design.mat, h, model,
        weights = weights, rawdata = rawdata, options = opt)
    if (!any(is.na(r$x))) {
        if (opt$band) {
            upper <- r$model.y + 2 * r$se
            upper <- pmin(pmax(upper, par()$usr[3]), par()$usr[4])
            lower <- r$model.y - 2 * r$se
            lower <- pmin(pmax(lower, par()$usr[3]), par()$usr[4])
            polygon(c(r$eval.points, rev(r$eval.points)), c(lower, rev(upper)),
                    col = 0, border = 0)
        }
        if (opt$display %in% "se") {
            upper <- r$estimate + 2 * r$se
            upper <- pmin(pmax(upper, par()$usr[3]), par()$usr[4])
            lower <- r$estimate - 2 * r$se
            lower <- pmin(pmax(lower, par()$usr[3]), par()$usr[4])
            lines(r$eval.points, upper, lty = 3, col = 0)
            lines(r$eval.points, lower, lty = 3, col = 0)
        }
        lines(r$eval.points, r$estimate, col = 0)
    }
    if (opt$band) {
        upper <- rnew$model.y + 2 * rnew$se
        upper <- pmin(pmax(upper, par()$usr[3]), par()$usr[4])
        lower <- rnew$model.y - 2 * rnew$se
        lower <- pmin(pmax(lower, par()$usr[3]), par()$usr[4])
        polygon(c(rnew$eval.points, rev(rnew$eval.points)), c(lower,
            rev(upper)), col = "cyan", border = 0)
    }
    lines(rnew$eval.points, rnew$estimate, lty = opt$lty, col = opt$col)
    if ((model == "none") & (opt$display %in% "se")) {
        upper <- rnew$estimate + 2 * rnew$se
        upper <- pmin(pmax(upper, par()$usr[3]), par()$usr[4])
        lower <- rnew$estimate - 2 * rnew$se
        lower <- pmin(pmax(lower, par()$usr[3]), par()$usr[4])
        lines(rnew$eval.points, upper, lty = 3, col = opt$col)
        lines(rnew$eval.points, lower, lty = 3, col = opt$col)
    }
    if (!opt$add)
        points(rawdata$x, rawdata$y, col = 1, pch = opt$pch,
               cex = 2/log(rawdata$nobs))
    box(col = 1, lty = 1)
    rnew
}
"plot2" <-
function (latitude2, longitude2, theta, phi)
{
    par(pch = "x")
    a <- (longitude2 * pi)/180
    b <- (latitude2 * pi)/180
    radtheta <- (theta * pi)/180
    radphi <- (phi * pi)/180
    xyzcheck <<- ((cos(a) * cos(b) * sin(radtheta) * cos(radphi)) +
        (sin(b) * sin(radphi)) - (sin(a) * cos(b) * cos(radphi) *
        cos(radtheta)))
    long2 <<- a[xyzcheck >= 0]
    lat2 <<- b[xyzcheck >= 0]
    if (length(lat2) == 0) {
        points(0, 0, type = "n")
        text(0.6, -1.2, labels = "Data set:")
        break
    }
    X <<- (cos(long2) * cos(lat2) * cos(radtheta)) + (sin(long2) *
        cos(lat2) * sin(radtheta))
    Y <<- (sin(lat2) * cos(radphi)) + ((cos(lat2) * sin(radphi)) *
        ((sin(long2) * cos(radtheta)) - (cos(long2) * sin(radtheta))))
    points(X, Y)
}
"plot2d" <-
function (d, f, theta, phi)
{
    par(pch = "*")
    a <- (f * pi)/180
    b <- (d * pi)/180
    radtheta <- (theta * pi)/180
    radphi <- (phi * pi)/180
    xyzcheck <<- ((cos(a) * cos(b) * sin(radtheta) * cos(radphi)) +
        (sin(b) * sin(radphi)) - (sin(a) * cos(b) * cos(radphi) *
        cos(radtheta)))
    llong <<- a[xyzcheck >= 0]
    llat <<- b[xyzcheck >= 0]
    invislong <<- a[xyzcheck < 0]
    invislat <<- b[xyzcheck < 0]
    if (length(llat) == 0) {
        par(pty = "s")
        plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "",
            xlim = c(-1, 1), ylim = c(-1, 1))
        list(invislong = invislong, invislat = invislat)
        break
    }
    X <<- (cos(llong) * cos(llat) * cos(radtheta)) + (sin(llong) *
        cos(llat) * sin(radtheta))
    Y <<- (sin(llat) * cos(radphi)) + ((cos(llat) * sin(radphi)) *
        ((sin(llong) * cos(radtheta)) - (cos(llong) * sin(radtheta))))
    par(pty = "s")
    plot(X, Y, axes = FALSE, xlab = "", ylab = "", xlim = c(-1, 1),
         ylim = c(-1, 1))
    list(invislong = invislong, invislat = invislat)
}
"print.graph" <-
function (file, ...)
{
    dev.print(file = file, ...)
    invisible()
}
"provide.data" <-
function (data, path, options = list())
{
    describe <- sm.options(options)$describe
    name <- deparse(substitute(data))
    if (missing(path))
        path <- system.file("smdata", package="sm")
    datafile <- file.path(path, paste(name, ".dat", sep = ""))
    docfile <- file.path(path, paste(name, ".doc", sep = ""))
    if (!exists(name, where=.GlobalEnv, inherits = FALSE)) {
        if (file.exists(datafile)) {
            cat("Data file being loaded\n")
            assign(name, read.table(datafile, header = TRUE),
                   envir = .GlobalEnv)
            attach(what = data, name = name)
        }
        else cat("Data file does not exist\n")
    }
    else {
        if (!is.data.frame(data))
            cat("object exists, not as a data.frame\n")
        else {
            cat(paste(name, "already loaded\n"))
            attach.frame(data, name = name)
        }
    }
    if (describe) {
        if(file.exists(docfile)) file.show(docfile)
        else cat("Doc file not found\n")
    }
    invisible()
}
"replace.na" <-
function (List, comp, value)
{
    arg <- paste(substitute(List), "$", substitute(comp), sep = "")
    arg.value <- eval(parse(text = arg), parent.frame(1))
    if (any(is.na(arg.value))) {
        change <- paste(arg, "<-", deparse(substitute(value)))
        a <- eval(parse(text = change), parent.frame(1))
    }
    invisible()
}
"sig.trace" <-
function (expn, hvec, ...)
{
    opt <- sm.options(list(...))
    replace.na(opt, display, "lines")
    expn.char <- paste(deparse(substitute(expn)), collapse = "")
    lead.char <- substring(expn.char, 1, nchar(expn.char) - 1)
    nvec <- length(hvec)
    pvec <- vector("numeric", length = nvec)
    for (i in 1:nvec) {
        extn.char <- paste(lead.char, ", h = ", as.character(hvec[i]), ")")
        result <- eval(parse(text = extn.char))
        pvec[i] <- result$p
    }
    if (!(opt$display %in% "none")) {
        plot(hvec, pvec, type = "l", ylim = c(0, max(pvec)),
            xlab = "Smoothing parameter, h", ylab = "p-value")
        if (max(pvec) >= 0.05)
            lines(range(hvec), c(0.05, 0.05), lty = 2)
    }
    invisible(list(h = hvec, p = pvec))
}
"sj" <-
function (x, h)
{
    phi6 <- function(x) (x^6 - 15 * x^4 + 45 * x^2 - 15) * dnorm(x)
    phi4 <- function(x) (x^4 - 6 * x^2 + 3) * dnorm(x)
    n <- length(x)
    lambda <- quantile(x, 0.75) - quantile(x, 0.25)
    a <- 0.92 * lambda * n^(-1/7)
    b <- 0.912 * lambda * n^(-1/9)
    W <- matrix(rep(x, rep(n, n)), ncol = n, byrow = TRUE)
    W <- W - matrix(rep(x, n), ncol = n, byrow = TRUE)
    W1 <- matrix(phi6(W/b), ncol = n)
    tdb <- as.numeric(rep(1, n) %*% W1 %*% rep(1, n))
    tdb <- -tdb/(n * (n - 1) * b^7)
    W1 <- matrix(phi4(W/a), ncol = n)
    sda <- as.numeric(rep(1, n) %*% W1 %*% rep(1, n))
    sda <- sda/(n * (n - 1) * a^5)
    alpha2 <- 1.357 * (abs(sda/tdb))^(1/7) * h^(5/7)
    W1 <- matrix(phi4(W/alpha2), ncol = n)
    sdalpha2 <- as.numeric(rep(1, n) %*% W1 %*% rep(1, n))
    sdalpha2 <- sdalpha2/(n * (n - 1) * alpha2^5)
    result <- (dnorm(0, sd = sqrt(2))/(n * abs(sdalpha2)))^0.2 - h
    attributes(result)$names <- NULL
    as.double(result)
}
"sm.ancova" <-
function (x, y, group, h, model = "none", band = TRUE, test = TRUE,
    h.alpha = 2 * diff(range(x))/length(x), weights = as.integer(rep(1,
        length(x))), ...)
{
    opt <- sm.options(list(...))
    x.name <- deparse(substitute(x))
    y.name <- deparse(substitute(y))
    if (any(is.na(c(x, y))) | any(is.na(group))) {
        xy <- cbind(x, y)
        ok <- as.logical(apply(!is.na(xy), 1, prod)) & (!is.na(group))
        xy <- xy[ok, ]
        y <- as.vector(xy[, ncol(xy)])
        x <- xy[, -ncol(xy), drop = TRUE]
        group <- group[ok]
        cat("warning: missing data are removed\n")
    }
    replace.na(opt, display, "lines")
    replace.na(opt, ngrid, 50)
    replace.na(opt, xlab, x.name)
    replace.na(opt, ylab, y.name)
    ndim <- 1
    nobs <- length(x)
    replace.na(opt, nbins, round((nobs > 500) * 8 * log(nobs)/ndim))
    if (!missing(weights)) {
        if (!is.na(opt$nbins) & opt$nbins != 0)
            stop("if weights are set, nbins must be 0 or NA")
        weights <- as.vector(weights)
        if (any(weights < 0 | is.na(weights)))
            stop("Negative or NA weights are meaningless")
        if (!isInteger(weights)) {
            weights <- round(weights/min(weights[weights > 0]))
            cat("Warning: weights have been rescaled to integer values\n")
        }
    }
    replace.na(opt, nbins, round((nobs > 500) * 8 * log(nobs)/ndim))
    if (model == "none") {
        band <- FALSE
        test <- FALSE
    }
    fact <- factor(group)
    ord <- order(fact, x)
    xx <- x[ord]
    yy <- y[ord]
    weights <- weights[ord]
    fact <- fact[ord]
    fac.levels <- levels(fact)
    nlevels <- length(fac.levels)
    rawdata <- list(x = xx, y = yy, fac = fact, nbins = opt$nbins,
        nobs = nobs, ndim = ndim, devs = 0)
    if (opt$nbins > 0) {
        for (i in 1:nlevels) {
            ind <- (fact == fac.levels[i])
            bins <- binning(xx[ind], yy[ind], nbins = opt$nbins)
            if (i == 1) {
                x <- bins$x
                y <- bins$means
                fac <- rep(fac.levels[1], length(bins$x))
                weights <- bins$x.freq
                rawdata$devs <- bins$devs
            }
            else {
                x <- c(x, bins$x)
                y <- c(y, bins$means)
                fac <- c(fac, rep(fac.levels[i], length(bins$x)))
                weights <- c(weights, bins$x.freq)
                rawdata$devs <- c(rawdata$devs, bins$devs)
            }
        }
        weights <- as.integer(weights)
        fac <- factor(fac)
    }
    else {
        x <- xx
        y <- yy
        fac <- fact
    }
    n <- table(fac)
    eval.points <- opt$eval.points
    if (any(is.na(eval.points))) {
        start.eval <- max(tapply(x, fac, min))
        stop.eval <- min(tapply(x, fac, max))
        eval.points <- seq(start.eval, stop.eval, length = opt$ngrid)
    }
    if (!(opt$display %in% "none")) {
        plot(rawdata$x, rawdata$y, type = "n", xlab = opt$xlab,
            ylab = opt$ylab)
        text(rawdata$x, rawdata$y, as.character(rawdata$fac))
        if (!opt$band) {
            for (i in 1:nlevels) {
                ind <- (fac == fac.levels[i])
                sm.regression(x[ind], y[ind], h = h, weights = weights[ind],
                  ngrid = opt$ngrid, add = TRUE, lty = i)
            }
        }
    }
    B <- diag(0, sum(n))
    Sd <- diag(0, sum(n))
    istart <- 1
    for (i in 1:nlevels) {
        irange <- istart:(istart + n[i] - 1)
        xi <- x[irange]
        wi <- weights[irange]
        B[irange, irange] <- sm.sigweight(xi, weights = wi, options = opt)
        Sd[irange, irange] <- sm.weight(xi, xi, h, weights = wi, options = opt)
        istart <- istart + n[i]
    }
    Ss <- sm.weight(x, x, h, weights = weights, options = opt)
    sigma <- sqrt((y %*% B %*% y)[1, 1] + sum(rawdata$devs))
    if (model == "equal") {
        Q <- Sd - Ss
        Q <- t(Q) %*% diag(weights) %*% Q
        obs <- ((y %*% Q %*% y)/sigma^2)[1, 1]
        covar <- diag(1/weights)
    }
    if (model == "parallel") {
        D <- matrix(0, ncol = nlevels - 1, nrow = sum(n))
        istart <- n[1] + 1
        for (i in 2:nlevels) {
            D[istart:(istart + n[i] - 1), i - 1] <- 1
        }
        Q <- diag(sum(n)) - sm.weight(x, x, h.alpha, weights = weights,
            options = opt)
        Q <- solve(t(D) %*% t(Q) %*% diag(weights) %*% Q %*%
            D) %*% t(D) %*% t(Q) %*% diag(weights) %*% Q
        alpha <- as.vector(Q %*% y)
        covar <- diag(1/weights)
        covar <- rbind(cbind(Q %*% covar %*% t(Q), Q %*% covar),
            cbind(t(Q %*% covar), covar))
        adjy <- y - D %*% alpha
        ghat <- Ss %*% adjy
        ghati <- Sd %*% y
        obs <- sum(weights * (as.vector(D %*% alpha) + ghat -
            ghati)^2)/sigma^2
        Q <- cbind((diag(sum(n)) - Ss) %*% D, (Ss - Sd))
        Q <- t(Q) %*% diag(weights) %*% Q
        B <- rbind(matrix(0, nrow = nlevels - 1, ncol = sum(n) +
            nlevels - 1), cbind(matrix(0, nrow = sum(n), ncol = nlevels -
            1), B))
    }
    p <- NULL
    if (!(model == "none")) {
        p <- p.quad.moment(Q - B * obs, covar, obs, sum(weights) -
            length(weights))
        cat("Test of ", model, " lines:   h = ", signif(h), "   p-value = ",
            round(p, 4), "\n")
    }
    sigma <- sigma/sqrt(sum(weights) - 2 * nlevels)
    if (opt$band & !(opt$display %in% "none")) {
        if (nlevels > 2)
            print("Band available only to compare two groups.")
        else {
            ind <- (fac == fac.levels[1])
            model1 <- sm.regression(x[ind], y[ind], h = h,
                                    eval.points = eval.points,
                weights = weights[ind], options = opt, display = "none",
                ngrid = opt$ngrid, add = TRUE, lty = 1)
            ind <- fac == fac.levels[2]
            model2 <- sm.regression(x[ind], y[ind], h = h,
                                    eval.points = eval.points,
                weights = weights[ind], options = opt, display = "none",
                ngrid = opt$ngrid, add = TRUE, lty = 2)
            model.y <- (model1$estimate + model2$estimate)/2
            if (model == "parallel") {
                model.y <- cbind(model.y - alpha/2, model.y +
                  alpha/2)
            }
            se <- sqrt((model1$se/model1$sigma)^2 + (model2$se/model2$sigma)^2)
            se <- se * sigma
            upper <- model.y + se
            lower <- model.y - se
            if (model == "equal") {
                upper <- pmin(pmax(upper, par()$usr[3]), par()$usr[4])
                lower <- pmin(pmax(lower, par()$usr[3]), par()$usr[4])
                polygon(c(eval.points, rev(eval.points)), c(lower,
                  rev(upper)), border = 0, col = "cyan")
            }
            else if (model == "parallel") {
                upper[, 1] <- pmin(pmax(upper[, 1], par()$usr[3]),
                  par()$usr[4])
                lower[, 1] <- pmin(pmax(lower[, 1], par()$usr[3]),
                  par()$usr[4])
                upper[, 2] <- pmin(pmax(upper[, 2], par()$usr[3]),
                  par()$usr[4])
                lower[, 2] <- pmin(pmax(lower[, 2], par()$usr[3]),
                  par()$usr[4])
                polygon(c(eval.points, rev(eval.points)), c(lower[,
                  1], rev(upper[, 1])), border = 0, col = 6)
                polygon(c(eval.points, rev(eval.points)), c(lower[,
                  2], rev(upper[, 2])), border = 0, col = 5)
            }
            text(rawdata$x, rawdata$y, as.character(rawdata$fac))
            lines(eval.points, model1$estimate, lty = 1)
            lines(eval.points, model2$estimate, lty = 2)
        }
    }
    r <- list(p = p, model = model, sigma = sigma)
    if (model == "parallel")
        r <- list(p = p, model = model, sigma = sigma, alphahat = alpha)
    r$data <- list(x = x, y = y, group = fac, nbins = rawdata$nbins,
        devs = rawdata$devs, weights = weights)
    r$call <- match.call()
    invisible(r)
}
"sm.autoregression" <-
function (x, h = hnorm(x), d = 1, maxlag = d, lags, se = FALSE, ask = TRUE)
{
    sm.autoregression.1d <- function(x, h, x.name, lags, se = FALSE,
        ask = FALSE) {
        n <- length(x)
        if (any(diff(lags)) < 0)
            stop("lags must be in increasing order")
        x2.name <- paste(x.name, "(t)", sep = "")
        xlow <- min(x) - diff(range(x))/20
        xhi <- max(x) + diff(range(x))/20
        lags <- sort(lags)
        for (m in lags) {
            x1 <- x[(m + 1):n]
            x0 <- x[1:(n - m)]
            r <- sm.regression.eval.1d(x0, x1, h = h, model = "none",
                options = list(hmult = 1))
            x1.name <- paste(x.name, "(t-", as.character(m),
                ")", sep = "")
            plot(x0, x1, xlim = c(xlow, xhi), ylim = c(xlow,
                xhi), xlab = x1.name, ylab = x2.name)
            lines(r$eval.points, r$estimate)
            if (se) {
                rho1 <- acf(x0, lag.max = 1, plot = FALSE)$acf[2]
                lines(r$eval.points, r$estimate + 2 * r$se/sqrt(1 - rho1),
                      lty = 3)
                lines(r$eval.points, r$estimate - 2 * r$se/sqrt(1 - rho1),
                      lty = 3)
            }
            title(paste("Regression of ", x.name, " on past data",
                        sep = ""))
            if (ask & (m < lags[length(lags)]))
                pause()
        }
        invisible(r)
    }
    sm.autoregression.2d <- function(x, h, x.name, lags, ask = ask,
        ngrid = 20, display = "none") {
        if (dim(lags)[2] != 2)
            stop("dim(lags)[2] must be 2")
        evpt <- seq(quantile(x, 0.1), quantile(x, 0.9), length = ngrid)
        n <- length(x)
        nplot <- dim(lags)[1]
        for (ip in 1:nplot) {
            m1 <- min(lags[ip, ])
            m2 <- max(lags[ip, ])
            x0 <- x[1:(n - m2)]
            x1 <- x[(m2 - m1 + 1):(n - m1)]
            x2 <- x[(m2 + 1):n]
            r <- sm.regression.eval.2d(cbind(x0, x1), x2, h = c(h,
                h), model = "none", eval.points = cbind(evpt,
                evpt), weights = rep(1, n - m2), options = list(hmult = 1,
                h.weights = rep(1, n - m2), poly.index = 1))
            persp(evpt, evpt, r)
            head <- paste("Regression of ", x.name, " on past data (lags: ",
                as.character(m1), ", ", as.character(m2), ")",
                sep = "")
            title(head)
            if (ask & (ip < nplot)) pause()
        }
        invisible(r)
    }
    x.name <- deparse(substitute(x))
    if (missing(lags)) {
        if (d == 1)
            lags <- (1:maxlag)
        else lags <- cbind(1:(maxlag - 1), 2:maxlag)
    }
    else {
        if (isMatrix(lags))
            d <- 2
    }
    x <- as.vector(x)
    if (d == 1)
        r <- sm.autoregression.1d(x, h, x.name, lags, se = se,
            ask = ask)
    else r <- sm.autoregression.2d(x, h, x.name, lags, ask = ask)
    invisible(r)
}
"sm.binomial" <-
function (x, y, N = rep(1, length(y)), h, ...)
{
    n <- length(y)
    x.name <- deparse(substitute(x))
    y.name <- deparse(substitute(y))
    opt <- sm.options(list(...))
    if (any(is.na(c(x, y)))) {
        xy <- cbind(x, y)
        ok <- as.logical(apply(!is.na(xy), 1, prod))
        xy <- xy[ok, ]
        y <- as.vector(xy[, ncol(xy)])
        x <- xy[, -ncol(xy), drop = TRUE]
        cat("warning: missing data are removed\n")
    }
    replace.na(opt, display, "estimate")
    replace.na(opt, ngrid, 25)
    replace.na(opt, ylim, c(0, 1))
    replace.na(opt, nbins, round((n > 100) * 8 * log(n)))
    display <- opt$display
    if (length(x) != n)
        stop("x and y have different length")
    if (length(N) != n)
        stop("N and y have different length")
    y <- as.integer(y)
    if (min(diff(x)) < 0) {
        y <- y[order(x)]
        x <- sort(x)
    }
    replace.na(opt, eval.points, seq(min(x), max(x), length = opt$ngrid))
    replace.na(opt, nbins, round((nobs > 100) * 8 * log(n)/ndim))
    if (all(N == 1))
        yplot <- jitter(y, amount = 0)
    else yplot <- y/N
    if (display != "none" & opt$add == FALSE) {
        replace.na(opt, xlab, x.name)
        replace.na(opt, ylab, paste("Pr{", y.name, "}", sep = ""))
        plot(x, yplot, ylim = opt$ylim, xlab = opt$xlab, ylab = opt$ylab,
            col = 1, type = "n")
        abline(0, 0, col = 1, lty = 3)
        abline(1, 0, col = 1, lty = 3)
    }
    if (display != "none")
        points(x, yplot, pch = opt$pch, col = opt$col)
    rawdata <- list(x = x, y = y, N = N, nbins = opt$nbins, nobs = n, ndim = 1)
    if (opt$nbins > 0) {
        bins <- binning(x, y, nbins = opt$nbins)
        binsN <- binning(x, N, nbins = opt$nbins)
        x <- bins$x
        y <- round(bins$sums)
        N <- round(binsN$sums)
        nx <- length(y)
    }
    result <- sm.glm(x, cbind(y, N - y), family = binomial(link = logit),
        h = h, eval.points = opt$eval.points, start = log((y +
            0.5)/(N - y + 1)))
    result$call <- match.call()
    if (display != "none") {
        lines(result$eval.points, result$estimate, col = opt$col)
        if (display == "se") {
            lines(result$eval.points, result$lower, lty = 3, col = opt$col)
            lines(result$eval.points, result$upper, lty = 3, col = opt$col)
        }
    }
    result$data <- list(x = x, y = y, N = N, nbins = opt$nbins)
    invisible(result)
}
"sm.binomial.bootstrap" <-
function (x, y, N = rep(1, length(x)), h, nboot = 99, degree = 1,
    fixed.disp = FALSE, family = binomial(logit), ...)
{
    rbetabinom <- function(n, size, prob, disp) {
        if (disp > 1 & min(size) > 2 & min(size) > disp) {
            psi <- (disp - 1)/(size - 1)
            alpha <- prob * (1/psi - 1)
            beta <- (1 - prob) * (1/psi - 1)
            p <- rbeta(n, alpha, beta)
            y <- rbinom(n, size, p)
        }
        else y <- rbinom(n, size, prob)
        return(y)
    }
    D <- function (mu, y, wt) sum(family$dev.resids(y, mu, wt))
    n <- length(x)
    sm <- sm.binomial(x, y, N, h, xlab = deparse(substitute(x)),
        ylab = paste("Pr{", deparse(substitute(y)), "}", sep = ""), ...)
    X <- cbind(1, poly(x, degree))
    colnames(X) <- seq(len=ncol(X))
    glm.model <- glm.fit(X, cbind(y, N - y), family = family, ...)
    glm.fitted <- fitted(glm.model)
    glm.resid <- residuals(glm.model)
    lines(x, glm.fitted, lty = 2, col = 2)
    p.boot <- 0
    sm.orig <- sm.binomial(x, y, N, h, eval.points = x, display = "none")
    sm.fitted <- sm.orig$estimate
    disp.orig <- D(sm.fitted, y/N, N)/(n - degree - 1)
    if (fixed.disp) disp <- 1
    else disp <- disp.orig
    ts.orig <- (D(glm.fitted, y/N, N) - D(sm.fitted, y/N, N))/disp
    cat("Dispersion parameter = ", disp.orig, "\n")
    cat("Test statistic = ", ts.orig, "\n")
    yboot <- rep(NA, n)
    for (i in 1:nboot) {
        yboot <- rbetabinom(n, N, glm.fitted, disp)
        cat("Sample:", i, " ")
        sm.fit <- sm.glm(x, cbind(yboot, N - yboot), family = family,
            h, eval.points = x, start = log((yboot + 0.5)/(N - yboot + 0.5)))
        sm.fitted <- sm.fit$estimate
        ts.boot <- (D(glm.fitted, yboot/N, N) - D(sm.fitted, yboot/N, N))/disp
        if (ts.boot > ts.orig) p.boot <- p.boot + 1
        lines(x, sm.fitted, lty = 2, col = 4)
    }
    cat("\n")
    lines(sm$eval.points, sm$estimate)
    p.boot <- p.boot/(nboot + 1)
    cat("Observed significance = ", p.boot, "\n")
    invisible(list(call = match.call(), significance = p.boot,
        test.statistic = ts.orig, dispersion = disp.orig))
}
"sm.density" <-
function (x, h, model = "none", weights = rep(1, nobs), ...)
{
    x.name <- deparse(substitute(x))
    opt <- sm.options(list(...))
    positive <- opt$positive
    if (any(is.na(x))) {
        ok <- as.logical(apply(!is.na(as.matrix(x)), 1, prod))
        if (is.vector(x)) x <- as.vector(x[ok])
        else x <- x[ok, ]
        cat("warning: missing data are removed\n")
    }
    if (length(dim(x)) > 0) {
        ndim <- dim(x)[2]
        if (ndim > 3) ndim <- 3
        nobs <- dim(x)[1]
    }
    else {
        ndim <- 1
        nobs <- length(x)
    }
    if (!missing(weights)) {
        if (!is.na(opt$nbins) & opt$nbins != 0)
            stop("if weights are set, nbins must be either 0 or NA")
        weights <- as.vector(weights)
        if (any(weights < 0 | is.na(weights)))
            stop("Negative or NA weights are meaningless")
        if (!isInteger(weights)) {
            weights <- round(weights/min(weights[weights > 0]))
            cat("Warning: weights have been rescaled to integer values\n")
        }
    }
    replace.na(opt, nbins, round((nobs > 500) * 8 * log(nobs)/ndim))
    if (nobs != length(weights))
        stop("length of x and weights do not match")
    X <- na.omit(cbind(x, weights))
    weights <- X[, ndim + 1]
    x <- X[, 1:ndim]
    if (ndim == 1) x <- as.vector(x)
    rawdata <- list(nbins = opt$nbins, x = x, nobs = nobs, ndim = ndim)
    if (opt$nbins > 0 & ndim < 3) {
        bins <- binning(x, nbins = opt$nbins)
        x <- bins$x
        weights <- bins$x.freq
        nx <- length(x)
        if (!all(is.na(opt$h.weights)))
            stop("use of h.weights is incompatible with binning - set nbins=0")
    }
    else nx <- nobs
    if (positive) {
        if (ndim == 1)
            replace.na(opt, delta, min(x))
        if (ndim == 2)
            replace.na(opt, delta, apply(x, 2, min))
    }
    if (missing(h)) {
        if (positive) {
            if (ndim == 1)
                h <- hnorm(log(x + opt$delta), weights)
            else if (ndim == 2)
                h <- hnorm(log(x + outer(rep(1, nx), opt$delta)), weights)
        }
        else h <- hnorm(x, weights)
    }
    if (ndim == 1) {
        if (length(h) != 1)
            stop("length(h) != 1")
        replace.na(opt, xlab, x.name)
        replace.na(opt, ylab, "Probability density function")
        replace.na(opt, delta, min(x))
        est <- sm.density.1d(x, h, model, weights, rawdata, options = opt)
    }
    if (ndim == 2) {
        if (length(h) != 2)
            stop("length(h) != 2")
        dimn <- dimnames(x)[[2]]
        name.comp <- if (!is.null(dimn) & !all(dimn == "")) dimn
        else {
            if (!is.null(attributes(x)$names))
                attributes(x)$names
            else outer(x.name, c("[1]", "[2]"), paste, sep = "")
        }
        replace.na(opt, xlab, name.comp[1])
        replace.na(opt, ylab, name.comp[2])
        replace.na(opt, delta, c(min(x[, 1]), min(x[, 2])))
        est <- sm.density.2d(x, h, weights, rawdata, options = opt)
    }
    if (ndim == 3) {
        dimn <- dimnames(x)[[2]]
        name.comp <- if (!is.null(dimn) & !all(dimn == "")) dimn
        else {
            if (!is.null(attributes(x)$names))
                attributes(x)$names
            else outer(x.name, c("[1]", "[2]", "[3]"), paste, sep = "")
        }
        replace.na(opt, xlab, name.comp[1])
        replace.na(opt, ylab, name.comp[2])
        replace.na(opt, zlab, name.comp[3])
        opt$nbins <- 0
        est <- sm.density.3d(x, h = h, contour = opt$props[1], options = opt)
    }
    est$data <- list(x = x, nbins = opt$nbins, freq = weights)
    est$call <- match.call()
    invisible(est)
}
"sm.density.1d" <-
function (x, h = hnorm(x, weights), model = "none", weights,
    rawdata = list(x = x), options = list())
{
    absent <- function(x) missing(x) | any(is.na(x) | is.null(x))
    opt <- sm.options(options)
    replace.na(opt, display, "estimate")
    panel <- opt$panel
    band <- opt$band
    hmult <- opt$hmult
    if (any(is.na(opt$h.weights)))
        replace.na(opt, h.weights, rep(1, length(x)))
    else band <- panel <- FALSE
    if (model == "none")
        band <- FALSE
    if (opt$add | opt$display %in% "none")
        panel <- FALSE
    a <- if (opt$positive) c(0, max(x) * 1.05)
    else c(min(x) - diff(range(x))/4, max(x) + diff(range(x))/4)
    replace.na(opt, xlim, a)
    long.x <- rep(x, weights)
    a <- if (opt$positive)
        max(0.4/(quantile(long.x, 0.75) - quantile(long.x, 0.5)),
            0.4/(quantile(long.x, 0.5) - quantile(long.x, 0.25)))
    else 0.6/sqrt(wvar(x, weights))
    replace.na(opt, yht, a)
    replace.na(opt, ylim, c(0, opt$yht))
    replace.na(opt, ngrid, 100)
    if (!opt$add & !(opt$display %in% "none"))
        plot(opt$xlim, opt$ylim, type = "n", xlab = opt$xlab, ylab = opt$ylab)
    opt$band <- band
    opt$panel <- panel
    if (!(opt$display %in% "none"))
        est <- smplot.density(x, h, weights, rawdata, options = opt)
    if (panel) {
        save.rug <- opt$rugplot
        items <- c("Bandwidth:", "  - normal optimal", "  - plug-in",
            "  - xval", "  - increase", "  - decrease", "  - movie up",
            "  - movie down", "Add Normal band", "Exit")
        if (rawdata$nbins > 0) items <- items[-c(3, 4)]
        hsj.flag <- FALSE
        hcv.flag <- FALSE
        ind <- menu(items, graphics = TRUE, title = "Density estimation")
        while (items[ind] != "Exit") {
            if (items[ind] == "  - normal optimal") {
                h <- hnorm(x, weights)
                hmult <- 1
            }
            if (items[ind] == "  - plug-in") {
                if (!hsj.flag) {
                  h.sj <- hsj(x)
                  hsj.flag <- TRUE
                }
                h <- h.sj
                hmult <- 1
            }
            if (items[ind] == "  - xval") {
                if (!hcv.flag) {
                  h.cv <- hcv(x)
                  hcv.flag <- TRUE
                }
                h <- h.cv
                hmult <- 1
            }
            else if (items[ind] == "  - increase") {
                hmult <- hmult * 1.1
            }
            else if (items[ind] == "  - decrease") {
                hmult <- hmult/1.1
            }
            else if (items[ind] == "  - movie up") {
                for (i in 1:6) {
                  hmult <- hmult * 1.1
                  opt$hmult <- hmult
                  opt$rugplot <- FALSE
                  est <- smplot.density(x, h, options = opt)
                }
                hmult <- hmult * 1.1
            }
            else if (items[ind] == "  - movie down") {
                for (i in 1:6) {
                  hmult <- hmult/1.1
                  opt$hmult <- hmult
                  opt$rugplot <- FALSE
                  est <- smplot.density(x, h, options = opt)
                }
                hmult <- hmult/1.1
            }
            else if (items[ind] == "Add Normal band" | items[ind] ==
                "Remove Normal band") {
                band <- !band
                if (items[ind] == "Add Normal band") {
                  items[ind] <- "Remove Normal band"
                }
                else (items[ind] <- "Add Normal band")
            }
            opt$hmult <- hmult
            opt$band <- band
            est <- smplot.density(x, h, options = opt)
            opt$rugplot <- save.rug
            cat("h = ", signif(h * hmult, 7), "\n")
            ind <- menu(items, graphics = TRUE, title = "Density estimation")
        }
    }
    if (all(!is.na(opt$eval.points))) {
        if (opt$positive)
            est <- sm.density.positive.1d(x, h, weights, options = opt)
        else est <- sm.density.eval.1d(x, h, weights, options = opt)
    }
    else if (opt$display %in% "none")
        est <- sm.density.eval.1d(x, h, weights = weights, options = opt)
    if (all(opt$h.weights == rep(1, length(x))) & opt$positive ==
        FALSE) {
        se <- sqrt(dnorm(0, sd = sqrt(2))/(4 * sum(weights) *
            h))
        upper <- sqrt(est$estimate) + 2 * se
        lower <- pmax(sqrt(est$estimate) - 2 * se, 0)
        upper <- upper^2
        lower <- lower^2
        est$se <- rep(se, length(est$eval.points))
        est$upper <- upper
        est$lower <- lower
    }
    invisible(est)
}
"sm.density.2d" <-
function (X, h = hnorm(X, weights), weights = rep(1, length(x)),
    rawdata = list(), options = list())
{
    x <- X[, 1]
    y <- X[, 2]
    opt <- sm.options(options)
    replace.na(opt, display, "persp")
    replace.na(opt, ngrid, 50)
    replace.na(opt, zlab, "Density function")
    replace.na(opt, xlim, range(X[, 1]))
    replace.na(opt, ylim, range(X[, 2]))
    display <- opt$display
    if (display == "none") opt$panel <- FALSE
    replace.na(opt, h.weights, rep(1, length(x)))
    hmult <- opt$hmult
    if (display == "persp")
        est <- sm.persplot(x, y, h, weights, rawdata, options = opt)
    else {
        if (display == "image")
            est <- sm.imageplot(x, y, h, weights, rawdata, options = opt)
        else if (display == "slice")
            est <- sm.sliceplot(x, y, h, weights, rawdata, options = opt)
    }
    if (opt$panel) {
        items <- c("Bandwidth:", "  - Normal optimal", "  - increase",
            "  - decrease", "Exit")
        ind <- menu(items, graphics = TRUE, title = "2-d density estimation")
        while (items[ind] != "Exit") {
            if (items[ind] == "  - Normal optimal") {
                hmult <- 1
            }
            else if (items[ind] == "  - increase") {
                hmult <- hmult * 1.1
            }
            else if (items[ind] == "  - decrease") {
                hmult <- hmult/1.1
            }
            if (display == "persp")
                est <- sm.persplot(x, y, h, weights, rawdata,
                  options = opt)
            else {
                if (display == "image")
                  est <- sm.imageplot(x, y, h, weights, rawdata, options = opt)
                else if (display == "slice")
                  est <- sm.sliceplot(x, y, h, weights, rawdata, options = opt)
            }
            ind <- menu(items, graphics = TRUE,
                        title = "2-d density estimation")
        }
    }
    if (isMatrix(opt$eval.points))
        est <- sm.density.eval.2d(x, y, h, xnew = opt$eval.points[,
            1], ynew = opt$eval.points[, 2], eval.type = "points",
            weights = weights, options = opt)
    else {
        if (display == "none")
            est <- sm.density.eval.2d(x, y, h, eval.type = "grid",
                weights = weights, options = opt)
    }
    if (all(opt$h.weights == rep(1, length(x)))) {
        se <- dnorm(0, sd = sqrt(2))/sqrt(4 * sum(weights) * h[1] * h[2])
        upper <- sqrt(est$estimate) + 2 * se
        lower <- pmax(sqrt(est$estimate) - 2 * se, 0)
        upper <- upper^2
        lower <- lower^2
        est$se <- est$estimate - est$estimate + se
        est$upper <- upper
        est$lower <- lower
    }
    invisible(est)
}
"sm.density.3d" <-
function (data, h = hnorm(data), contour = 75, shadow = TRUE, colour = TRUE,
    title.colour = 3, label.colour = 1, axes.colour = 1, plot.colour = 1,
    shadow.colour = 2, plt = TRUE, cex = NA, maxpoly = 10000, options = list())
{
    opt <- sm.options(options)
    replace.na(opt, ngrid, 20)
    replace.na(opt, hmult, 1)
    hmult <- opt$hmult
    xlab <- opt$xlab
    ylab <- opt$ylab
    zlab <- opt$zlab
    theta <- opt$theta
    phi <- opt$theta
    ngrid <- opt$ngrid
    x <- data[, 1]
    y <- data[, 2]
    z <- data[, 3]
    if (is.na(xlab))
        if (!is.null(attributes(data)$dimnames)) {
            xlab <- attributes(data)$dimnames[[2]][1]
            ylab <- attributes(data)$dimnames[[2]][2]
            zlab <- attributes(data)$dimnames[[2]][3]
        }
        else {
            xlab <- ""
            ylab <- ""
            zlab <- ""
        }
    lng <- nrow(data)
    xmin <- min(x)
    xmax <- max(x)
    ymin <- min(y)
    ymax <- max(y)
    zmin <- min(z)
    zmax <- max(z)
    xbar <- mean(x)
    ybar <- mean(y)
    zbar <- mean(z)
    if (is.na(ngrid))
        ngrid <- 12
    gridsize <- ngrid
    hx <- h[1] * hmult
    hy <- h[2] * hmult
    hz <- h[3] * hmult
    fhat <- vector()
    for (i in (1:lng)) {
        fhat[i] <- sum(exp(-0.5 * (((x[i] - x)/hx)^2 + ((y[i] -
            y)/hy)^2 + ((z[i] - z)/hz)^2)))
    }
    fhat <- fhat/(lng * hx * hy * hz * (2 * pi)^(3/2))
    fhat <- sort(fhat)
    fmax <- fhat[lng]
    p <- (lng * (100 - contour))/100
    p.trunc <- trunc(p)
    height <- fhat[p.trunc] + (p - p.trunc) *
        (fhat[p.trunc + 1] - fhat[p.trunc])
    height <- (height * 90)/fmax
    xx <- seq(xmin, xmax, length = gridsize)
    yy <- seq(ymin, ymax, length = gridsize)
    zz <- seq(zmin, zmax, length = gridsize)
    tmp <- as.double(rep(0, maxpoly * 3))
    result <- .Fortran("npcont", as.double(x), as.double(y),
        as.double(z), as.double(xx), as.double(yy), as.double(zz),
        as.double(hx), as.double(hy), as.double(hz), as.integer(lng),
        as.integer(gridsize), as.double(fmax), as.double(height),
        as.integer(maxpoly), as.double(rep(0, gridsize^3)), as.integer(0),
        as.double(rep(0, maxpoly * 3 * 3)), tmp, tmp, tmp,
                       PACKAGE="sm")
    lng.coord <- 3 * result[[16]]
    xcoord <- result[[18]][1:lng.coord]
    ycoord <- result[[19]][1:lng.coord]
    zcoord <- result[[20]][1:lng.coord]
    coord <- cbind(xcoord, ycoord, zcoord)
    attributes(coord)$dimnames <- list(character(), c(xlab, ylab,
        zlab))
    if (plt) {
        np.contour.plot.3d.(coord, data, shadow, gridsize, 3,
            xmin, xmax, ymin, ymax, zmin, zmax, xlab, ylab, zlab,
            theta, phi, colour, title.colour, label.colour, axes.colour,
            plot.colour, shadow.colour, cex)
    }
    else coord
}
"sm.density.compare" <-
function (x, group, h = NA, model = "none", test = TRUE, nboot = 100,
    monitor = TRUE, ...)
{
    opt <- sm.options(list(...))
    fact <- factor(group)
    fact.levels <- levels(fact)
    nlev <- length(fact.levels)
    ni <- table(fact)

    replace.na(opt, ngrid, 50)
    replace.na(opt, display, "lines")
    replace.na(opt, xlab, deparse(substitute(x)))
    replace.na(opt, ylab, "Density")
    replace.na(opt, xlim, c(min(x) - diff(range(x))/4, max(x) +
                            diff(range(x))/4))
    replace.na(opt, eval.points,
               seq(opt$xlim[1], opt$xlim[2], length=opt$ngrid))
    if(length(opt$lty < nlev)) opt$lty <- 1:nlev

    band <- opt$band
    ngrid <- opt$ngrid
    xlim <- opt$xlim
    y <- x
    if (model == "none") {
        band <- FALSE
        test <- FALSE
    }
    if (opt$display %in% "none")
        band <- FALSE
    if (band & (nlev > 2)) {
        cat("Reference band available to compare two groups only.", "\n")
        band <- FALSE
    }
    if (is.na(h)) {
        hvec <- tapply(y, fact, FUN = "hnorm")
        h <- exp(mean(log(hvec)))
    }
    opt$band <- band
    opt$test <- test
    estimate <- matrix(0, ncol = opt$ngrid, nrow = nlev)
    se <- matrix(0, ncol = opt$ngrid, nrow = nlev)
    for (i in 1:nlev) {
        sm <- sm.density(y[fact == fact.levels[i]], h = h, display = "none",
                         eval.points = opt$eval.points)
        estimate[i, ] <- sm$estimate
        se[i, ] <- sm$se
    }
    eval.points <- sm$eval.points
    if (!(opt$display %in% "none" | band)) {
        plot(xlim, c(0, 1.1 * max(as.vector(estimate))), xlab = opt$xlab,
            ylab = opt$ylab, type = "n")
        for (i in 1:nlev) lines(eval.points, estimate[i, ],
            lty = opt$lty[i])
    }
    est <- NULL
    p <- NULL
    if (model == "equal" & test) {
        if (nlev == 2) {
            ts <- sum((estimate[1, ] - estimate[2, ])^2)
        }
        else {
            sm.mean <- sm.density(y, h = h, xlim = opt$xlim,
                ngrid = opt$ngrid, display = "none")$estimate
            ts <- 0
            for (i in 1:nlev) ts <- ts + ni[i] *
                sum((estimate[i,] - sm.mean)^2)
        }
        p <- 0
        est.star <- matrix(0, ncol = opt$ngrid, nrow = nlev)
        for (i in 1:nboot) {
            ind <- (1:length(y))
            for (i in 1:nlev) {
                indi <- sample((1:length(ind)), ni[i])
                est.star[i, ] <- sm.density(y[ind[indi]], h = h,
                  ngrid = opt$ngrid, xlim = opt$xlim, display = "none")$estimate
                ind <- ind[-indi]
            }
            if (nlev == 2) {
                ts.star <- sum((est.star[1, ] - est.star[2, ])^2)
            }
            else {
                sm.mean <- sm.density(y, h = h, xlim = opt$xlim,
                  ngrid = opt$ngrid, display = "none")$estimate
                ts.star <- 0
                for (i in 1:nlev) {
                  ts.star <- ts.star + ni[i] * sum((est.star[i,] - sm.mean)^2)
                }
            }
            if (ts.star > ts)
                p <- p + 1
            if (monitor) {
                cat(i)
                cat(" ")
            }
        }
        p <- p/nboot
        cat("\nTest of equal densities:  p-value = ", round(p,3), "\n")
        est <- list(p = p, h = h)
    }
    if (model == "equal" & band) {
        av <- (sqrt(estimate[1, ]) + sqrt(estimate[2, ]))/2
        se <- sqrt(se[1, ]^2 + se[2, ]^2)
        upper <- (av + se)^2
        lower <- pmax(av - se, 0)^2
        plot(xlim, c(0, 1.1 * max(as.vector(estimate), upper)),
            xlab = opt$xlab, ylab = opt$ylab, type = "n")
        polygon(c(eval.points, rev(eval.points)), c(upper, rev(lower)),
            col = "cyan", border = 0)
        lines(eval.points, estimate[1, ], lty = opt$lty[1])
        lines(eval.points, estimate[2, ], lty = opt$lty[2])
        est <- list(p = p, upper = upper, lower = lower, h = h)
    }
    invisible(est)
}
"sm.density.eval.1d" <-
function (x, h, weights = rep(1, n), options = list())
{
    opt <- sm.options(options)
    replace.na(opt, h.weights, rep(1, length(x)))
    replace.na(opt, xlim, c(min(x) - diff(range(x))/4, max(x) +
        diff(range(x))/4))
    replace.na(opt, ngrid, 100)
    hmult <- opt$hmult
    h.weights <- opt$h.weights
    xlim <- opt$xlim
    ngrid <- opt$ngrid
    replace.na(opt, eval.points, seq(xlim[1], xlim[2], length = ngrid))
    xnew <- opt$eval.points
    n <- length(x)
    neval <- length(xnew)
    W <- matrix(rep(xnew, rep(n, neval)), ncol = n, byrow = TRUE)
    W <- W - matrix(rep(x, neval), ncol = n, byrow = TRUE)
    W1 <- matrix(rep(h.weights, neval), ncol = n, byrow = TRUE)
    W <- exp(-0.5 * (W/(hmult * h * W1))^2)/W1
    est <- W %*% weights/(sum(weights) * sqrt(2 * pi) * hmult * h)
    invisible(list(eval.points = xnew, estimate = as.vector(est),
        h = h * hmult, h.weights = h.weights, weights = weights))
}
"sm.density.eval.2d" <-
function (x, y, h, xnew, ynew, eval.type = "points", weights = rep(1, n),
          options = list())
{
    opt <- sm.options(options)
    replace.na(opt, xlim, range(x))
    replace.na(opt, ylim, range(y))
    replace.na(opt, ngrid, 50)
    replace.na(opt, h.weights, rep(1, length(x)))
    if (missing(xnew))
        xnew <- seq(opt$xlim[1], opt$xlim[2], length = opt$ngrid)
    if (missing(ynew))
        ynew <- seq(opt$ylim[1], opt$ylim[2], length = opt$ngrid)
    n <- length(x)
    nnew <- length(xnew)
    h.weights <- opt$h.weights
    hmult <- opt$hmult
    W1 <- matrix(rep(xnew, rep(n, nnew)), ncol = n, byrow = TRUE)
    W1 <- W1 - matrix(rep(x, nnew), ncol = n, byrow = TRUE)
    W2 <- matrix(rep(h.weights, nnew), ncol = n, byrow = TRUE)
    Wx <- exp(-0.5 * (W1/(hmult * h[1] * W2))^2)/W2
    W1 <- matrix(rep(ynew, rep(n, nnew)), ncol = n, byrow = TRUE)
    W1 <- W1 - matrix(rep(y, nnew), ncol = n, byrow = TRUE)
    Wy <- exp(-0.5 * (W1/(opt$hmult * h[2] * W2))^2)/W2
    if (eval.type == "points")
        est <- as.vector(((Wx * Wy) %*% weights)/(sum(weights) *
            2 * pi * h[1] * h[2] * hmult^2))
    else est <- (Wx %*% (weights * t(Wy)))/(sum(weights) * 2 *
        pi * h[1] * h[2] * hmult^2)
    invisible(list(eval.points = cbind(xnew, ynew), estimate = est,
        h = h * hmult, h.weights = h.weights, weights = weights))
}
"sm.density.positive.1d" <-
function (x, h, weights, options = list())
{
    if (min(x) <= 0)
        cat("Warning: some data are not positive\n")
    opt <- sm.options(options)
    delta <- opt$delta
    replace.na(opt, ngrid, 100)
    replace.na(opt, xlim, c(0, max(x)))
    if (min(opt$xlim) < 0)
        cat("Warning: xlim<0 with positive=TRUE \n")
    if (missing(h))
        h <- hnorm(log(x + delta), weights)
    ngrid <- opt$ngrid
    ev.pt <- opt$eval.points
    if (any(is.na(ev.pt))) {
        a <- log(opt$xlim + 1/ngrid)
        ev.pt <- exp(seq(min(a), max(a), length = ngrid)) - 1/ngrid
    }
    opt$eval.points <- log(ev.pt + delta)
    f <- sm.density.eval.1d(log(x + delta), h = h, weights = weights,
        options = opt)
    est <- f$estimate/(ev.pt + delta)
    est[is.na(est)] <- 0
    list(eval.points = ev.pt, estimate = as.vector(est), h = h)
}
"sm.density.positive.2d" <-
function (X, h = c(hnorm(log(X[, 1] + delta[1]), weights), hnorm(log(X[,
    2] + delta[2]), weights)), eval.type = "points", weights = rep(1, nrow(X)),
          options = list())
{
    opt <- sm.options(options)
    replace.na(opt, ngrid, 50)
    replace.na(opt, delta, apply(X, 2, min))
    if (min(X) <= 0)
        cat("Warning: some data are not positive\n")
    if (dim(X)[2] != 2)
        cat("parameter X must be two-columns matrix\n")
    x1 <- X[, 1]
    x2 <- X[, 2]
    delta <- opt$delta
    replace.na(opt, xlim, range(x1))
    replace.na(opt, ylim, range(x2))
    replace.na(opt, ngrid, 50)
    xlim <- opt$xlim
    ylim <- opt$ylim
    ngrid <- opt$ngrid
    ax <- log(xlim + 1/ngrid)
    ay <- log(ylim + 1/ngrid)
    eval1 <- exp(seq(ax[1], ax[2], length = ngrid)) - 1/ngrid
    eval2 <- exp(seq(ay[1], ay[2], length = ngrid)) - 1/ngrid
    replace.na(opt, eval.points, cbind(eval1, eval2))
    eval1 <- opt$eval.points[, 1]
    eval2 <- opt$eval.points[, 2]
    pdf <- sm.density.eval.2d(log(x1 + delta[1]), log(x2 + delta[2]),
        h = h, xnew = log(eval1 + delta[1]), ynew = log(eval2 +
            delta[2]), eval.type = eval.type, weights = weights)
    if (eval.type == "points")
        est <- pdf$estimate/((eval1 + delta[1]) * (eval2 + delta[2]))
    else est <- pdf$estimate/outer(eval1 + delta[1], eval2 +
        delta[2])
    invisible(list(x1 = eval1, x2 = eval2, estimate = est, h = h))
}
"sm.density.positive.grid" <-
function (X, h = c(hnorm(log(X[, 1] + delta[1])), hnorm(log(X[,
    2] + delta[2]))), delta = c(min(X[, 1]), min(X[, 2])), ngrid = 50,
    eval.points = NA, xlim = range(X[, 1]), ylim = range(X[,
        2]))
{
    f <- sm.density.positive.2d(X, h, delta, ngrid = ngrid, xlim = xlim,
        ylim = ylim, eval.type = "grid")
    xx <- rep(f$x1, length(f$x2))
    yy <- rep(f$x2, rep(length(f$x2), length(f$x1)))
    zz <- as.vector(f$est, byrow = TRUE)
    f.int <- interp(xx, yy, zz)
    invisible(list(eval.points = cbind(f.int$x, f.int$y), estimate = f.int$z,
        h = h))
}
"sm.glm" <-
function (x, y, family, h, eval.points, start, offset)
{
    n <- length(x)
    X <- cbind(rep(1, n + 1), c(x, 0))
    ## in R, avoid zero weight
    if (isMatrix(y)) Y <- rbind(y, rep(1, ncol(y)))
    else Y <- c(y, 0)
    start <- c(start, 0)
    neval <- length(eval.points)
    if (missing(offset)) offset <- rep(0, n)
    W <- matrix(rep(eval.points, rep(n, neval)), ncol = n, byrow = TRUE)
    W <- W - matrix(rep(X[1:n, 2], neval), ncol = n, byrow = TRUE)
    W <- exp(-0.5 * (W/h)^2)
    cat("Cycles per point: ")
    est <- LP <- st.err <- dev <- var.eta <- rep(NA, neval)
    for (k in 1:neval) {
        X[n + 1, 2] <- eval.points[k]
        colnames(X) <- 1:2
        ## weight 0 does not work in R
        fit <- glm.fit(X, Y, w = c(W[k, ], 1e-8), family = family,
            etastart = start, offset = c(offset, 0))
        start <- fit$linear.predictors
        LP[k] <- start[n + 1]
        dev[k] <- fit$deviance
        cat(fit$iter)
        cat(" ")
        s <- W[k, ]
#         if (family$family[1] == "Binomial")
#             w <- apply(Y, 1, sum)[1:n]
#         else w <- rep(1, n)
        mu <- fit$fitted.values[1:n]
        Wk <- diag(s * fit$weights[1:n])
        XXinv <- solve(t(X[1:n, ]) %*% Wk %*% X[1:n, ])
        Li <- XXinv %*% t(X[1:n, ]) %*% diag(s)
        var.Bi <- Li %*% diag(fit$weights[1:n]) %*% t(Li)
        var.eta[k] <- t(X[n + 1, ]) %*% var.Bi %*% as.vector(X[n + 1, ])
    }
    cat("\n")
    st.err <- sqrt(var.eta)
    est <- family$linkinv(LP)
    result <- list(call = match.call(), eval.points = eval.points,
        estimate = est, lower = family$linkinv(LP - 2 * st.err),
        upper = family$linkinv(LP + 2 * st.err), linear.predictor = LP,
        se = st.err, deviance = dev)
    invisible(result)
}
"sm.imageplot" <-
function (x, y, h, weights, rawdata, options = list())
{
    opt <- sm.options(options)
    replace.na(opt, h.weights, rep(1, length(x)))
    replace.na(opt, xlab, deparse(substitute(x)))
    replace.na(opt, ylab, deparse(substitute(y)))
    replace.na(opt, zlab, "Density function")
    replace.na(opt, ngrid, 50)
    replace.na(opt, xlim, range(x))
    replace.na(opt, ylim, range(y))
    ngrid <- opt$ngrid
    xlim <- opt$xlim
    ylim <- opt$ylim
    xgrid <- seq(xlim[1], xlim[2], length = ngrid)
    ygrid <- seq(ylim[1], ylim[2], length = ngrid)
    if (!opt$positive)
        dgrid <- sm.density.eval.2d(x, y, h, xgrid, ygrid, eval.type = "grid",
            weights, opt)$estimate
    else {
        f <- sm.density.positive.grid(cbind(x, y), h, opt$delta,
            NA, opt$ngrid, opt$xlim, opt$ylim)
        xgrid <- f$eval.points[, 1]
        ygrid <- f$eval.points[, 2]
        dgrid <- f$estimate
    }
    image(xgrid, ygrid, dgrid, xlab = opt$xlab, ylab = opt$ylab)
    invisible(list(eval.points = cbind(xgrid, ygrid), estimate = dgrid,
        h = h * opt$hmult, h.weights = opt$h.weights, weights = weights))
}

"sm.options" <-
function (...)
{
    if (nargs() == 0) return(.sm.Options)
    current <- .sm.Options
    temp <- list(...)
    if (length(temp) == 1 && is.null(names(temp))) {
        arg <- temp[[1]]
        switch(mode(arg),
               list = temp <- arg,
               character = return(.sm.Options[arg]),
               stop(paste("invalid argument:", arg)))
    }
    if (length(temp) == 0) return(current)
    n <- names(temp)
    if (is.null(n)) stop("options must be given by name")
    changed <- current[n]
    current[n] <- temp
    if (sys.parent() == 0) env <- asNamespace("sm") else env <- parent.frame()
    assign(".sm.Options", current, envir = env)
    invisible(current)
}

"sm.persplot" <-
function (x, y, h = hnorm(cbind(x, y), weights), weights, rawdata = list(),
    options = opt)
{
    opt <- sm.options(options)
    replace.na(opt, h.weights, rep(1, length(x)))
    replace.na(opt, xlab, deparse(substitute(x)))
    replace.na(opt, ylab, deparse(substitute(y)))
    replace.na(opt, zlab, "Density function")
    replace.na(opt, ngrid, 50)
    replace.na(opt, xlim, range(x))
    replace.na(opt, ylim, range(y))
    ngrid <- opt$ngrid
    xlim <- opt$xlim
    ylim <- opt$ylim
    xgrid <- seq(xlim[1], xlim[2], length = ngrid)
    ygrid <- seq(ylim[1], ylim[2], length = ngrid)
    if (!opt$positive)
        dgrid <- sm.density.eval.2d(x, y, h, xgrid, ygrid, eval.type = "grid",
            weights, options = opt)$estimate
    else {
        f <- sm.density.positive.grid(cbind(x, y), h, opt$delta,
            ngrid, NA, xlim, ylim)
        xgrid <- f$eval.points[, 1]
        ygrid <- f$eval.points[, 2]
        dgrid <- f$estimate
    }
    persp(xgrid, ygrid, dgrid, xlab = opt$xlab, ylab = opt$ylab,
        zlab = opt$zlab, theta = -30, phi = 40, d = 4)
    invisible(list(eval.points = cbind(xgrid, ygrid), estimate = dgrid,
        h = h * opt$hmult, h.weights = opt$h.weights, weights = weights))
}
"sm.poisson" <-
function (x, y, h, ...)
{
    x.name <- deparse(substitute(x))
    y.name <- deparse(substitute(y))
    opt <- sm.options(list(...))
    if (any(is.na(c(x, y)))) {
        xy <- cbind(x, y)
        ok <- as.logical(apply(!is.na(xy), 1, prod))
        xy <- xy[ok, ]
        y <- as.vector(xy[, ncol(xy)])
        x <- xy[, -ncol(xy), drop = TRUE]
        cat("warning: missing data are removed\n")
    }
    family <- poisson(link = log)
    y <- as.integer(y)
    n <- length(y)
    replace.na(opt, display, "estimate")
    replace.na(opt, ngrid, 25)
    replace.na(opt, ylim, c(0, 1))
    replace.na(opt, pch, 1)
    replace.na(opt, col, 2)
    replace.na(opt, nbins, round((n > 100) * 8 * log(n)))
    display <- opt$display
    if (min(diff(x)) < 0) {
        y <- y[order(x)]
        x <- sort(x)
    }
    if (!(opt$display %in% "none") & opt$add %in% FALSE) {
        replace.na(opt, xlab, x.name)
        replace.na(opt, ylab, y.name)
        plot(x, y, xlab = opt$xlab, ylab = opt$ylab, col = 1, type = "n")
    }
    if (display != "none")
        points(x, y, pch = opt$pch, col = opt$col)
    replace.na(opt, eval.points, seq(min(x), max(x), length = opt$ngrid))
    rawdata <- list(x = x, y = y, nbins = opt$nbins, nobs = n,
        ndim = 1)
    if (opt$nbins > 0) {
        bins <- binning(x, y, nbins = opt$nbins)
        x <- bins$x
        y <- round(bins$sums)
        nx <- length(y)
        freq <- bins$x.freq
    }
    else freq <- rep(1, n)
    result <- sm.glm(x, y, family = family, h = h, eval.points = opt$eval.points,
        start = log(pmax(0.167, y)), offset = log(freq))
    result$call <- match.call()
    if (display != "none") {
        lines(result$eval.points, result$estimate, col = opt$col)
        if (display == "se") {
            lines(result$eval.points, result$lower, lty = 3, col = opt$col)
            lines(result$eval.points, result$upper, lty = 3, col = opt$col)
        }
    }
    result$data <- list(x = x, y = y, weights = freq, nbins = opt$nbins)
    invisible(result)
}
"sm.poisson.bootstrap" <-
function (x, y, h, nboot = 99, degree = 1, fixed.disp = FALSE,
          intercept = TRUE, family = poisson(link = log), ...)
{
    rNegBin <- function(n, mean, disp) {
        if (disp > 1) {
            p <- 1/disp
            r <- mean/(disp - 1)
            theta <- (rgamma(n, r) * (1 - p))/p
            y <- rpois(n, theta)
        }
        else y <- rpois(n, mean)
        return(y)
    }
    D <- function(mu, y, w, residuals = FALSE)
        sum(family$dev.resids(y, mu, w))

    x.name <- deparse(substitute(x))
    y.name <- deparse(substitute(y))
    y <- as.integer(y)
    y <- y[order(x)]
    x <- sort(x)
    sm <- sm.poisson(x, y, h, xlab = x.name, ylab = y.name, col = 3, ...)
    if (intercept)
        X <- cbind(1, poly(x, degree))
    else X <- outer(x, 1:degree, "^")
    colnames(X) <- seq(len=ncol(X))
    glm.model <- glm.fit(X, y, family = family)
    glm.fitted <- fitted(glm.model)
    lines(x, glm.fitted, col = 5)
    p.boot <- 0
    sm.orig <- sm.poisson(x, y, h, eval.points = x, display = "none",
        ...)
    sm.fitted <- sm.orig$estimate
    disp.orig <- D(sm.fitted, y, 1)/(length(y) - ncol(X))
    if (fixed.disp)
        disp <- 1
    else disp <- disp.orig
    ts.orig <- (D(glm.fitted, y, 1) - D(sm.fitted, y, 1))/disp
    cat("Dipersion parameter = ", disp.orig, "\n")
    cat("Test statistic = ", ts.orig, "\n")
    for (i in 1:nboot) {
        cat(i, " ")
        yboot <- rNegBin(length(glm.fitted), glm.fitted, disp)
        sm <- sm.poisson(x, yboot, h, eval.points = x, display = "none")
        sm.fitted <- sm$estimate
        ts.boot <- (D(glm.fitted, yboot, 1) - D(sm.fitted, yboot, 1))/disp
        if (ts.boot > ts.orig)
            p.boot <- p.boot + 1
        lines(x, sm.fitted, lty = 2, col = 6)
    }
    lines(sm$eval.points, sm$estimate, col = 3)
    lines(x, glm.fitted, col = 5)
    p.boot <- p.boot/(nboot + 1)
    cat("Observed significance = ")
    cat(p.boot)
    cat("\n")
    invisible(list(call = match.call(), test.statistic = ts.orig,
        significance = p.boot, disp = disp.orig))
}
"sm.regression" <-
function (x, y, h, design.mat = NA, model = "none", test = TRUE,
    weights = rep(1, nobs), ...)
{
    x.name <- deparse(substitute(x))
    y.name <- deparse(substitute(y))
    opt <- sm.options(list(...))
    if (any(is.na(c(x, y)))) {
        xy <- cbind(x, y)
        ok <- as.logical(apply(!is.na(xy), 1, prod))
        xy <- xy[ok, ]
        y <- as.vector(xy[, ncol(xy)])
        x <- xy[, -ncol(xy), drop = TRUE]
        cat("warning: missing data are removed\n")
    }
    if (length(dim(x)) > 0) {
        ndim <- 2
        nobs <- dim(x)[1]
    }
    else {
        ndim <- 1
        x <- as.vector(x)
        nobs <- length(x)
    }
    if (length(x)%/%ndim != length(y))
        stop("size of x and y do not match")
    if (!missing(h) & length(h) != ndim)
        stop("length(h) does not match size of x")
    replace.na(opt, display, "lines")
    replace.na(opt, poly.index, 1)
    replace.na(opt, nbins, round((nobs > 500) * 8 * log(nobs)/ndim))
    if (!missing(weights)) {
        if (!is.na(opt$nbins) & opt$nbins != 0)
            stop("if weights are set, nbins must be 0 or NA")
        weights <- as.vector(weights)
        if (any(weights < 0 | is.na(weights)))
            stop("Negative or NA weights are meaningless")
        if (!isInteger(weights)) {
            weights <- round(weights/min(weights[weights > 0]))
            cat("Warning: weights have been rescaled to integer values\n")
        }
    }
    rawdata <- list(x = x, y = y, nbins = opt$nbins, nobs = nobs, ndim = ndim)
    if (opt$nbins > 0) {
        bins <- binning(x, y, nbins = opt$nbins)
        x <- bins$x
        y <- bins$means
        weights <- bins$x.freq
        rawdata$devs <- bins$devs
        nx <- length(y)
    }
    else nx <- nobs
    replace.na(opt, h.weights, rep(1, nx))
    if (ndim == 1) {
        replace.na(opt, xlab, x.name)
        replace.na(opt, ylab, y.name)
        replace.na(opt, ngrid, 50)
        est <- sm.regression.1d(x, y, h, design.mat, model, weights,
            rawdata, options = opt)
    }
    else {
        replace.na(opt, ngrid, 20)
        dimn <- dimnames(x)[[2]]
        name.comp <- if (!is.null(dimn) & !all(dimn == "")) dimn
        else {
            if (!is.null(attributes(x)$names))
                attributes(x)$names
            else outer(x.name, c("[1]", "[2]"), paste, sep = "")
        }
        replace.na(opt, xlab, name.comp[1])
        replace.na(opt, ylab, name.comp[2])
        replace.na(opt, zlab, y.name)
        est <- sm.regression.2d(x, y, h, model, weights, rawdata,
                                options = opt)
    }
    est$data <- list(x = x, y = y, opt$nbins, freq = weights)
    est$call <- match.call()
    invisible(est)
}
"sm.regression.1d" <-
function (x, y, h, design.mat = NA, model = "none", weights = rep(1,
    length(x)), rawdata, options = list())
{
    opt <- sm.options(options)
    replace.na(opt, ngrid, 50)
    replace.na(opt, display, "lines")
    hmult <- opt$hmult
    if (model == "none") {
        opt$band <- FALSE
        opt$test <- FALSE
    }
    if (opt$add | opt$display %in% "none")
        opt$panel <- FALSE
    if (!(model == "none") & opt$panel == FALSE)
        opt$test <- TRUE
    r <- list(x = NA, y = NA, model.y = NA, se = NA, sigma = NA,
        h = h * hmult, hweights = opt$h.weights, weights = weights)
    if (!opt$add & !(opt$display %in% "none"))
        plot(rawdata$x, rawdata$y, xlab = opt$xlab, ylab = opt$ylab,
            type = "n")
    if (!(opt$display %in% "none")) {
        opt1 <- opt
        opt1$test <- FALSE
        r <- plot.regression(x, y, design.mat, h, r, model, weights,
            rawdata, options = opt1)
    }
    if (opt$test)
        rtest <- sm.regression.test(x, y, design.mat, h, model,
            weights, rawdata, options = opt)
    if (opt$panel) {
        items <- c("Bandwidth:", "  - increase", "  - decrease",
            "  - movie up", "  - movie down", "Exit")
        ind <- menu(items, graphics = TRUE, title = "Nonparametric regression")
        while (items[ind] != "Exit") {
            if (items[ind] == "  - increase") {
                hmult <- hmult * 1.1
            }
            else if (items[ind] == "  - decrease") {
                hmult <- hmult/1.1
            }
            else if (items[ind] == "  - movie up") {
                for (i in 1:6) {
                  hmult <- hmult * 1.1
                  opt1 <- opt
                  opt1$hmult <- hmult
                  opt1$test <- FALSE
                  r <- plot.regression(x, y, design.mat, h, r,
                    model, weights = weights, rawdata, options = opt1)
                }
                opt1$hmult <- hmult <- hmult * 1.1
            }
            else if (items[ind] == "  - movie down") {
                for (i in 1:6) {
                  hmult <- hmult/1.1
                  opt1 <- opt
                  opt1$hmult <- hmult
                  opt1$test <- FALSE
                  r <- plot.regression(x, y, design.mat, h, r,
                    model, weights = weights, rawdata, options = opt1)
                }
                opt1$hmult <- hmult <- hmult/1.1
            }
            else if (items[ind] == "Add linear band" | items[ind] ==
                "Remove linear band") {
                bandflag <- !bandflag
                if (!bandflag)
                  polygon(c(r$x, rev(r$x)), c(mean(y) - 2 * r$se,
                    mean(y) + 2 * rev(r$se)), col = 0)
                if (items[ind] == "Add linear band") {
                  items[ind] <- "Remove linear band"
                }
                else (items[ind] <- "Add linear band")
            }
            opt1 <- opt
            opt1$test <- FALSE
            opt1$hmult <- hmult
            r <- plot.regression(x, y, design.mat, h, r, model,
                weights = weights, rawdata, options = opt1)
            cat("h = ", signif(h * hmult, 7), "\n")
            ind <- menu(items, graphics = TRUE,
                        title = "Nonparametric regression")
        }
    }
    if (!(any(is.na(opt$eval.points))))
        r <- sm.regression.eval.1d(x, y, design.mat, h, model,
            weights, rawdata, options = opt)
    else if ((opt$display %in% "none") & (model == "none")) {
        opt$eval.points <- seq(min(x), max(x), length = opt$ngrid)
        r <- sm.regression.eval.1d(x, y, design.mat, h, model,
            weights, rawdata, options = opt)
    }
    if (opt$test)
        r <- list(eval.points = r$eval.points, estimate = r$estimate,
            model.y = r$model.y, se = r$se, sigma = r$sigma,
            h = r$h, hweights = r$hweights, weights = weights,
            model = rtest$model, p = rtest$p)
    r
}
"sm.regression.2d" <-
function (x, y, h, model = "none", weights = rep(1, n), rawdata,
    options = list())
{
    opt <- sm.options(options)
    replace.na(opt, h.weights, rep(1, length(y)))
    replace.na(opt, display, "persp")
    replace.na(opt, ngrid, 20)
    if (model == "none")
        opt$test <- FALSE
    n <- length(y)
    ev.points <- NA
    est <- NA
    x1grid <- seq(min(x[, 1]), max(x[, 1]), length = opt$ngrid)
    x2grid <- seq(min(x[, 2]), max(x[, 2]), length = opt$ngrid)
    e.points <- cbind(x1grid, x2grid)
    if (!(opt$display %in% "none")) {
        est <- sm.regression.eval.2d(x, y, h, model, e.points,
            opt$hull, weights, options = opt)
        ev.points <- e.points
        eye <- c(x1grid[1] + opt$eye.mult[1] * diff(range(x1grid)),
            x2grid[1] + opt$eye.mult[2] * diff(range(x2grid)),
            max(est[!is.na(est)]) + opt$eye.mult[3] * diff(range(est[!is.na(est)])))
        persp(x1grid, x2grid, est, xlab = opt$xlab, ylab = opt$ylab,
            zlab = opt$zlab, theta = -30, phi = 40, d = 4)
    }
    if (!(any(is.na(as.vector(opt$eval.points))))) {
        eval.type <- "points"
        ev.points <- opt$eval.points
        w <- sm.weight2(x, opt$eval.points, h, weights = weights,
            options = opt)
        est <- as.vector(w %*% y)
    }
    else if ((opt$display %in% "none") & (model == "none")) {
        est <- sm.regression.eval.2d(x, y, h, model, e.points,
            opt$hull, weights, options = opt)
        ev.points <- e.points
    }
    model.y <- NA
    se <- NA
    sigma <- NA
    r <- list(eval.points = ev.points, estimate = est, model.y = model.y,
        se = se, sigma = sigma, h = h * opt$hmult, hweights = opt$h.weights,
        weights = weights)
    if (opt$test) {
        rtest <- sm.regression.test(x, y, design.mat = NA, h,
            model, weights, rawdata, opt)
        r <- list(eval.points = ev.points, estimate = est, model.y = model.y,
            se = se, sigma = sigma, h = h * opt$hmult, hweights = opt$h.weights,
            model = rtest$model, p = rtest$p)
    }
    r
}
"sm.regression.autocor" <-
function (x = 1:n, y, h.first, minh, maxh, method = "direct",
    ...)
{
    GCV <- function(h, x, y, R, sqrt.R) {
        W <- sm.weight(x, x, h, options = list(hmult = 1))
        r <- (y - W %*% as.matrix(y))
        rss <- sum(r^2)
        Trace <- sum(diag(W))
        gcv.0 <- rss/(1 - Trace/length(x))^2
        Trace <- sum(diag(W %*% R))
        gcv.r <- rss/(1 - Trace/length(x))^2
        rw <- backsolve(sqrt.R, r)
        Trace <- sum(diag(W))
        gcv.ri <- sum(rw^2)/(1 - Trace/length(x))^2
        c(gcv.0, gcv.r, gcv.ri)
    }
    opt <- sm.options(list(...))
    replace.na(opt, display, "plot")
    replace.na(opt, ngrid, 15)
    ngrid <- opt$ngrid
    n <- length(y)
    if (length(x) != n)
        stop("x and y must have equal length\n")
    if (missing(minh) & missing(x))
        minh <- 0.5
    if (missing(maxh) & missing(x))
        maxh <- 10
    w <- sm.weight(x, x, h = h.first, options = list(hmult = 1))
    ym <- as.vector(w %*% y)
    r <- (y - ym)
    autocov <- rep(0, n)
    for (k in 0:2) {
        u <- r[1:(n - k)] * r[(k + 1):n]
        autocov[k + 1] <- sum(u)/n
    }
    var <- autocov[1]
    rho1 <- autocov[2]/var
    rho2 <- autocov[3]/var
    a1 <- rho1 * (1 - rho2)/(1 - rho1^2)
    a2 <- (rho2 - rho1^2)/(1 - rho1^2)
    cat("AR[1:2] coeff: ", c(a1, a2), "\n")
    for (k in 3:(n - 1)) autocov[k + 1] <- a1 * autocov[k] +
        a2 * autocov[k - 1]
    autocorr <- autocov/var
    R <- diag(n)
    R <- outer(1:n, 1:n, function(i, j, r) r[abs(i - j) + 1],
        r = autocorr)
    sqrt.R <- chol(R)
    hvector <- seq(minh, maxh, length = ngrid)
    min.gcv <- Inf
    h.opt <- 0
    result <- matrix(0, ngrid, 3, dimnames = list(NULL, c("no.cor",
        "direct", "indirect")))
    cat(paste("Search for h (runs up to ", as.character(ngrid),
        "): ", sep = "", collapse = NULL))
    for (i in 1:ngrid) {
        h <- hvector[i]
        result[i, ] <- GCV(h, x, y, R, sqrt.R)
        cat(" ")
        cat(i)
    }
    cat("\n")
    if (!(opt$display %in% "none")) {
        maxlag <- min(30, n - 1)
        acf <- array(autocorr[1:(maxlag + 1)], dim = c(maxlag +
            1, 1, 1))
        lag <- array(0:maxlag, dim = c(maxlag + 1, 1, 1))
#        acf.plot(list(acf = acf, lag = lag, type = "correlation",
#            series = "residuals from preliminary smoothing", n.used = n))
        plot(lag, acf, sub="residuals from preliminary smoothing", type="h")
        pause()
        plot(c(hvector[1], hvector[ngrid]), c(min(result), max(result)),
            type = "n", xlab = "h", ylab = "Generalised cross-validation")
        title(paste("GCV criterion, method:", method, collapse = NULL))
        lines(hvector, result[, method], col = 2)
        pause()
    }
    h1 <- hvector[order(result[, method])[1]]
    cat("Suggested value of h: ", h1, "\n")
    sm1 <- sm.regression.eval.1d(x, y, h = h1, model = "none",
        options = list(hmult = 1))
    if (missing(x))
        x.name <- "time"
    else x.name <- deparse(substitute(x))
    if (!(opt$display %in% "none")) {
        plot(x, y, xlab = x.name, ylab = deparse(substitute(y)),
            ...)
        lines(sm1$eval.points, sm1$estimate, col = 2)
    }
    sm1$aux <- list(h.first = h.first, first.sm = ym, acf = autocorr,
        raw.residuals = r)
    invisible(sm1)
}
"sm.regression.eval.1d" <-
function (x, y, design.mat, h, model = "none", weights = rep(1,
    length(x)), rawdata, options = list())
{
    opt <- sm.options(options)
    replace.na(opt, band, FALSE)
    replace.na(opt, test, FALSE)
    replace.na(opt, ngrid, 50)
    replace.na(opt, eval.points, seq(min(x), max(x), length = opt$ngrid))
    if (missing(rawdata))
        rawdata <- list(x = x, y = y, nbins = 0)
    band <- opt$band
    test <- opt$test
    ngrid <- opt$ngrid
    h.weights <- opt$h.weights
    eval.points <- opt$eval.points
    w <- sm.weight(x, eval.points, h, weights = weights, options = opt)
    est <- as.vector(w %*% y)
    sig <- sm.sigma(x, y, rawdata = rawdata, weights = weights)
    n <- length(x)
    ne <- length(eval.points)
    if (model == "none") {
        model.y <- est
        se <- as.vector(sig * sqrt(((w^2) %*% (1/weights))))
    }
    else if ((model == "no.effect") | (model == "no effect")) {
        if (is.na(as.vector(design.mat)[1])) {
            X <- matrix(rep(1, n), ncol = 1)
            model.y <- rep(wmean(y, weights), ne)
        }
        else {
            X <- design.mat
            model.y <- rep(0, ne)
        }
        X <- diag(n) - X %*% solve(t(X) %*% diag(weights) %*%
            X) %*% t(X) %*% diag(weights)
        se <- sig * sqrt(diag(w %*% X %*% diag(1/weights) %*%
            t(w)))
    }
    else if (model == "linear") {
        e <- cbind(rep(1, ne), eval.points - mean(x))
        l <- cbind(rep(1, n), x - mean(x))
        l <- e %*% solve(t(l) %*% diag(weights) %*% l) %*% t(l) %*%
            diag(weights)
        model.y <- as.vector(l %*% y)
        se <- as.vector(sig * sqrt(((w - l)^2) %*% (1/weights)))
    }
    list(eval.points = eval.points, estimate = est, model.y = model.y,
        se = se, sigma = sig, h = h * opt$hmult, hweights = h.weights,
        weights = weights)
}
"sm.regression.eval.2d" <-
function (x, y, h, model, eval.points, hull = TRUE, weights, options = list())
{
    opt <- sm.options(options)
    hmult <- opt$hmult
    h.weights <- opt$h.weights
    n <- nrow(x)
    ngrid <- nrow(eval.points)
    wd1 <- matrix(rep(eval.points[, 1], n), ncol = n)
    wd1 <- wd1 - matrix(rep(x[, 1], ngrid), ncol = n, byrow = TRUE)
    wd2 <- matrix(rep(eval.points[, 2], n), ncol = n)
    wd2 <- wd2 - matrix(rep(x[, 2], ngrid), ncol = n, byrow = TRUE)
    wy <- matrix(rep(h.weights, ngrid), ncol = n, byrow = TRUE)
    w1 <- exp(-0.5 * (wd1/(h[1] * hmult * wy))^2)
    w1 <- w1 * matrix(rep(weights, ngrid), ncol = n, byrow = TRUE)
    w2 <- exp(-0.5 * (wd2/(h[2] * hmult * wy))^2)
    wy <- matrix(rep(y, ngrid), ncol = n, byrow = TRUE)
    if (opt$poly.index == 0)
        est <- w1 %*% t(w2 * wy)/(w1 %*% t(w2))
    if (opt$poly.index == 1) {
        a11 <- w1 %*% t(w2)
        a12 <- (w1 * wd1) %*% t(w2)
        a13 <- w1 %*% t(w2 * wd2)
        a22 <- (w1 * wd1^2) %*% t(w2)
        a23 <- (w1 * wd1) %*% t(w2 * wd2)
        a33 <- w1 %*% t(w2 * wd2^2)
        d <- a22 * a33 - a23^2
        b1 <- 1/(a11 - ((a12 * a33 - a13 * a23) * a12 + (a13 *
            a22 - a12 * a23) * a13)/d)
        b2 <- (a13 * a23 - a12 * a33) * b1/d
        b3 <- (a12 * a23 - a13 * a22) * b1/d
        c1 <- w1 %*% t(w2 * wy)
        c2 <- (w1 * wd1) %*% t(w2 * wy)
        c3 <- w1 %*% t(w2 * wy * wd2)
        est <- b1 * c1 + b2 * c2 + b3 * c3
    }
    if (hull) {
        hull.points <- x[order(x[, 1], x[, 2]), ]
        dh <- diff(hull.points)
        hull.points <- hull.points[c(TRUE, !((dh[, 1] == 0) & (dh[,
            2] == 0))), ]
        hull.points <- hull.points[chull(hull.points), ]
        nh <- nrow(hull.points)
        gstep <- matrix(rep(eval.points[2, ] - eval.points[1,
            ], nh), ncol = 2, byrow = TRUE)
        hp.start <- matrix(rep(eval.points[1, ], nh), ncol = 2,
            byrow = TRUE)
        hull.points <- hp.start + gstep * round((hull.points -
            hp.start)/gstep)
        hull.points <- hull.points[chull(hull.points), ]
        grid.points <- cbind(rep(eval.points[, 1], ngrid), rep(eval.points[,
            2], rep(ngrid, ngrid)))
        D <- diff(rbind(hull.points, hull.points[1, ]))
        temp <- D[, 1]
        D[, 1] <- D[, 2]
        D[, 2] <- (-temp)
        C <- as.vector((hull.points * D) %*% rep(1, 2))
        C <- matrix(rep(C, ngrid^2), nrow = ngrid^2, byrow = TRUE)
        D <- t(D)
        wy <- ((grid.points %*% D) >= C)
        wy <- apply(wy, 1, all)
        wy[wy] <- 1
        wy[!wy] <- NA
        wy <- matrix(wy, ncol = ngrid)
    }
    else {
        w1 <- (w1 > exp(-2))
        w2 <- (w2 > exp(-2))
        wy <- w1 %*% t(w2)
        wy[wy > 0] <- 1
        wy[wy == 0] <- NA
    }
    est <- est * wy
    invisible(est)
}
"sm.regression.test" <-
function (x, y, design.mat = NA, h, model = "no.effect", weights = rep(1,
    length(y)), rawdata, options = list())
{
    opt <- sm.options(options)
    if (length(dim(x)) > 0) {
        ndim <- 2
        n <- dim(x)[1]
        W <- sm.weight2(x, x, h, weights = weights, option = opt)
        S <- cbind(rep(1, n), x[, 1] - mean(x[, 1]), x[, 2] -
            mean(x[, 2]))
    }
    else {
        ndim <- 1
        n <- length(x)
        W <- sm.weight(x, x, h, weights = weights, options = opt)
        S <- cbind(rep(1, n), x - mean(x))
    }
    if ((model == "no.effect") | (model == "no effect")) {
        if (is.na(as.vector(design.mat)[1]))
            S <- matrix(rep(1, n), ncol = 1)
        else S <- design.mat
    }
    if ((model == "linear") | (model == "no.effect") | (model ==
        "no effect")) {
        S <- diag(n) - S %*% solve(t(S) %*% diag(weights) %*%
            S) %*% t(S) %*% diag(weights)
        W <- diag(n) - W
        W <- t(W) %*% diag(weights) %*% W
        e <- as.vector(S %*% y)
        r0 <- sum(weights * e^2) + sum(rawdata$devs)
        r1 <- as.numeric(t(e) %*% W %*% e) + sum(rawdata$devs)
        ts <- (r0 - r1)/r1
        p <- p.quad.moment(diag(weights) - (1 + ts) * W, S %*%
            diag(1/weights), ts, sum(weights) - length(weights))
    }
    print(paste("Test of", model, "model:  significance = ",
        round(p, 3)))
    list(model = model, p = p, h = h * opt$hmult, hweights = opt$h.weights)
}
"sm.rm" <-
function (Time, y, minh = 0.1, maxh = 2, optimize = FALSE,
          rice.display = FALSE, ...)
{
    rice <- function(h, nSubj, Time, ym, var, r, poly.index = 1) {
        nTime <- length(Time)
        w <- sm.weight(Time, Time, h, options = list(poly.index = poly.index))
        fitted <- w %*% ym
        rss <- sum((ym - fitted)^2)
        Trace <- sum(diag(w %*% r))
        criterion <- sqrt(rss/nTime - (var/nSubj) * (1 - 2 *
            Trace/nTime))
        criterion
    }
    if (!isMatrix(y))
        stop("y must be a matrix")
    opt <- sm.options(list(...))
    replace.na(opt, ngrid, 20)
    ngrid <- opt$ngrid
    nSubj <- dim(y)[1]
    nTime <- dim(y)[2]
    if (missing(Time)) Time <- 1:nTime
    ym <- apply(y, 2, mean)
    z <- y - matrix(ym, nrow = nSubj, ncol = nTime, byrow = TRUE)
    autocov <- rep(0, nTime)
    for (k in 0:(nTime - 1)) {
        u <- z[, 1:(nTime - k)] * z[, (k + 1):nTime]
        autocov[k + 1] <- sum(u)/(nSubj * nTime)
    }
    var <- autocov[1]
    autocorr <- autocov/var
    cat("Autocovariances & autocorrelations:\n")
    print(matrix(cbind(autocov, autocorr), ncol = 2, dimnames = list(0:(nTime -
        1), c("auto-cov", "auto-corr"))))
    r <- diag(nTime)
    for (k in 1:nTime) {
        for (j in 1:nTime) r[k, j] <- autocorr[abs(k - j) + 1]
    }
    hvector <- seq(minh, maxh, length = ngrid)
    min.obj <- Inf
    h.opt <- 0
    cat("       Rice's criterion:\n")
    cat("       h    indept.   depend.\n")
    result <- matrix(0, ngrid, 2, dimnames = list(NULL, c("indept",
        "depend")))
    for (i in 1:ngrid) {
        h <- hvector[i]
        obj.0 <- rice(h, nSubj, Time, ym, var, diag(nTime), opt$poly.index)
        obj.r <- rice(h, nSubj, Time, ym, var, r, opt$poly.index)
        result[i, 1] <- obj.0
        result[i, 2] <- obj.r
        if (obj.r < min.obj) {
            min.obj <- obj.r
            h.opt <- h
        }
        print(c(h, obj.0, obj.r))
    }
    if (rice.display) {
        plot(c(hvector[1], hvector[ngrid]), c(min(result), max(result)),
            type = "n", xlab = "h", ylab = "sqrt(rice criterion)")
        title(main = "Modified Rice criterion for selecting h",
            sub = paste("dashed line: assume independence,",
                " continuous: allow for correlation", collapse = NULL))
        lines(hvector, result[, 1], lty = 3)
        lines(hvector, result[, 2], lty = 1)
        pause()
    }
    if (optimize) {
        cat("Search for optimum h using optim...\n")
        optimum <- optim(par = h.opt, fn = rice, method = "L-BFGS",
            lower = 0, nSubj = nSubj, Time = Time, ym = ym, var = var, r = r)
        print(optimum$par)
        h.opt <- optimum$par
    }
    cat("h: ", h.opt, "\n")
    if (opt$display %in% "se")
        display1 <- "lines"
    else display1 <- opt$display
    sm <- sm.regression(Time, ym, h = h.opt, hmult = 1, display = display1,
        ylab = paste(deparse(substitute(y)), "(mean values)",
            collapse = NULL), add = opt$add)
    if (opt$display %in% "se") {
        W <- sm.weight(Time, sm$eval.points, h = h.opt, options = list())
        V <- (var/nSubj) * r
        se <- sqrt(diag(W %*% V %*% t(W)))
        lines(sm$eval.points, sm$estimate + 2 * se, lty = 3)
        lines(sm$eval.points, sm$estimate - 2 * se, lty = 3)
    }
    sm$aux <- list(mean = ym, var = var, autocorr = autocorr, h = h.opt)
    invisible(sm)
}
"sm.sigma" <-
function (x, y, rawdata, weights = rep(1, length(y)), diff.ord = 2)
{
    n <- length(x)
    if (diff.ord == 1) {
        yd <- diff(y[order(x)])
        ww <- 1/weights[order(x)]
        wd <- ww[2:n] + ww[1:(n - 1)]
        ssq1 <- sum(yd^2/wd)
        ssq2 <- sum(rawdata$devs)
        sig <- sqrt((ssq1 + ssq2)/(sum(weights) - 1))
    }
    else {
        yy <- y[order(x)]
        xx <- sort(x)
        xx1 <- diff(xx)
        xx2 <- diff(xx, lag = 2)
        a <- xx1[-1]/xx2
        b <- xx1[-(n - 1)]/xx2
        a[xx2 == 0] <- 0.5
        b[xx2 == 0] <- 0.5
        ww <- weights[order(x)]
        cc <- a^2/ww[1:(n - 2)] + b^2/ww[3:n] + 1/ww[2:(n - 1)]
        eps <- yy[1:(n - 2)] * a + yy[3:n] * b - yy[2:(n - 1)]
        ssq1 <- sum(eps^2/cc)
        ssq2 <- sum(rawdata$devs)
        sig <- sqrt((ssq1 + ssq2)/(sum(ww) - 2))
    }
    sig
}
"sm.sigweight" <-
function (x, weights = rep(1, length(x)), ...)
{
    n <- length(x)
    xx <- sort(x)
    xx1 <- diff(xx)
    xx2 <- diff(xx, lag = 2)
    a <- xx1[-1]/xx2
    b <- xx1[-(n - 1)]/xx2
    a[xx2 == 0] <- 0.5
    b[xx2 == 0] <- 0.5
    c <- sqrt(a^2/weights[1:(n - 2)] + b^2/weights[3:n] + 1/weights[2:(n -
        1)])
    D <- cbind(rep(0, n - 2), diag(-1/c), rep(0, n - 2)) + cbind(diag(a/c),
        rep(0, n - 2), rep(0, n - 2)) + cbind(rep(0, n - 2),
        rep(0, n - 2), diag(b/c))
    D <- rbind(rep(0, n), D, rep(0, n))
    t(D) %*% D
}
"sm.sliceplot" <-
function (x, y, h, weights, rawdata = list(), options = list())
{
    opt <- sm.options(options)
    if (opt$positive) {
        cat("sliceplot not available with option positive=TRUE\n")
        cat("choose display='image' or display='persp'\n")
        stop()
    }
    replace.na(opt, h.weights, rep(1, length(x)))
    replace.na(opt, xlab, deparse(substitute(x)))
    replace.na(opt, ylab, deparse(substitute(y)))
    replace.na(opt, zlab, "Density function")
    replace.na(opt, ngrid, 50)
    replace.na(opt, xlim, range(x))
    replace.na(opt, ylim, range(y))
    ngrid <- opt$ngrid
    xlim <- opt$xlim
    ylim <- opt$ylim
    xgrid <- seq(xlim[1], xlim[2], length = ngrid)
    ygrid <- seq(ylim[1], ylim[2], length = ngrid)
    if (!opt$add) {
        plot(x, y, xlim = opt$xlim, ylim = opt$ylim, xlab = opt$xlab,
            ylab = opt$ylab, type = "n")
        points(rawdata$x[, 1], rawdata$x[, 2], col = 1,
            pch = opt$pch, cex = 2/log(rawdata$nobs))
    }
    if (opt$positive)
        f <- sm.density.positive.grid(cbind(x, y), h, opt$delta,
            NA, opt$ngrid, opt$xlim, opt$ylim)
    else f <- sm.density.eval.2d(x, y, h, xgrid, ygrid, eval.type = "grid",
        weights = weights)
    dgrid <- f$estimate
    xgrid <- f$eval.points[, 1]
    ygrid <- f$eval.points[, 2]
    if (opt$positive) {
        opt$eval.points <- cbind(x, y)
        dobs <- sm.density.positive.2d(cbind(x, y), h, eval.type = "points",
            weights = weights, options = opt)$estimate
    }
    else dobs <- sm.density.eval.2d(x, y, h, xnew = x, ynew = y,
        weights = weights)$estimate
    props <- opt$props
    hts <- quantile(rep(dobs, weights), prob = (100 - props)/100)
    for (i in 1:length(props)) {
        scale <- props[i]/hts[i]
        contour(xgrid, ygrid, dgrid * scale, level = hts[i] *
            scale, add = TRUE, lty = opt$lty, col = opt$col)
    }
    invisible(list(eval.points = cbind(xgrid, ygrid), estimate = dgrid,
        h = h * opt$hmult, h.weights = opt$h.weights, weights = weights))
}
"sm.sphere" <-
function (lat, long, kappa = 20, hidden = FALSE, sphim = FALSE,
          addpoints = FALSE, ...)
{
    opt <- sm.options(list(...))
    replace.na(opt, ngrid, 32)
    ngrid <- opt$ngrid
    panel <- opt$panel
    phi <- opt$phi
    theta <- opt$theta
    kap <- kappa
    invis <- plot2d(lat, long, theta, phi)
    sphdraw(theta, phi)
    if (!opt$panel) {
        if (hidden)
            hidplot(invis, theta, phi)
        if (sphim)
            sphimage(lat, long, kap, theta, phi, ngrid)
        if (sphim & addpoints)
            addplot(lat, long, theta, phi)
    }
    else {
        items <- c("Set theta and phi", "  - increase theta",
            "  - decrease theta", "  - increase phi", "  - decrease phi",
            "Add hidden points", "Add density estimate", "  - increase s.p.",
            "  - decrease s.p.", "  - add data points", "Exit")
        ind <- menu(items, graphics = TRUE, title = "Sphere")
        while (items[ind] != "Exit") {
            if (items[ind] == "Set theta and phi") {
                a <- change(theta, phi)
                theta <- a$theta
                phi <- a$phi
                invis <- plot2d(lat, long, theta, phi)
                sphdraw(theta, phi)
            }
            else if (items[ind] == "  - increase theta") {
                theta <- inctheta(theta, 30)
                invis <- plot2d(lat, long, theta, phi)
                sphdraw(theta, phi)
            }
            else if (items[ind] == "  - decrease theta") {
                theta <- inctheta(theta, -30)
                invis <- plot2d(lat, long, theta, phi)
                sphdraw(theta, phi)
            }
            else if (items[ind] == "  - increase phi") {
                phi <- incphi(phi, 30)
                invis <- plot2d(lat, long, theta, phi)
                sphdraw(theta, phi)
            }
            else if (items[ind] == "  - decrease phi") {
                phi <- incphi(phi, -30)
                invis <- plot2d(lat, long, theta, phi)
                sphdraw(theta, phi)
            }
            else if (items[ind] == "Add hidden points") {
                hidplot(invis, theta, phi)
            }
            else if (items[ind] == "Add density estimate") {
                sphimage(lat, long, kap, theta, phi, ngrid)
            }
            else if (items[ind] == "  - increase s.p.") {
                kap <- kap * 2
                sphimage(lat, long, kap, theta, phi, ngrid)
            }
            else if (items[ind] == "  - decrease s.p.") {
                kap <- kap/2
                sphimage(lat, long, kap, theta, phi, ngrid)
            }
            else if (items[ind] == "  - add data points") {
                par(pch = "*")
                addplot(lat, long, theta, phi)
            }
            else if (items[ind] == "Add 2nd data set") {
                par(pch = "x")
                addplot(lat2, long2, theta, phi)
            }
            ind <- menu(items, graphics = TRUE, title = "Sphere")
        }
    }
    par(pty = "m")
    invisible(list(theta = theta, phi = phi, kappa = kap))
}
"sm.survival" <-
function (x, y, status, h, hv = 0.05, p = 0.5, status.code = 1,
    ...)
{
    opt <- sm.options(list(...))
    replace.na(opt, display, "lines")
    replace.na(opt, ngrid, 50)
    replace.na(opt, xlab, deparse(substitute(x)))
    replace.na(opt, ylab, deparse(substitute(y)))
    replace.na(opt, eval.points, seq(min(x), max(x), length = opt$ngrid))
    eval.points <- opt$eval.points
    if (!(opt$display %in% "none" | opt$add == TRUE)) {
        plot(x, y, type = "n", xlab = opt$xlab, ylab = opt$ylab, ...)
        text(x[status == status.code], y[status == status.code], "x")
        text(x[status != status.code], y[status != status.code], "o")
    }
    n <- length(x)
    ne <- length(eval.points)
    xr <- x[order(y)]
    statusr <- status[order(y)]
    yr <- sort(y)
    w <- matrix(rep(eval.points, rep(n, ne)), ncol = n, byrow = TRUE)
    w <- w - matrix(rep(xr, ne), ncol = n, byrow = TRUE)
    w <- exp(-0.5 * (w/h)^2)
    wf <- t(apply(w, 1, rev))
    wf <- t(apply(wf, 1, cumsum))
    wf <- t(apply(wf, 1, rev))
    w <- w/wf
    st <- rep(0, n)
    st[statusr == status.code] <- 1
    w <- 1 - w * matrix(rep(st, ne), ncol = n, byrow = TRUE)
    w <- w[, st == 1]
    if (ne == 1)
        w <- matrix(w, ncol = length(w))
    yw <- yr[st == 1]
    w <- t(apply(w, 1, cumprod))
    w <- cbind(rep(1, ne), w)
    j <- -t(apply(w, 1, diff))
    J <- t(apply(j, 1, cumsum))
    wd <- J - p
    w <- exp(-0.5 * (wd/hv)^2)
    ns <- length(yw)
    s0 <- w %*% rep(1, ns)
    s1 <- (w * wd) %*% rep(1, ns)
    s2 <- (w * wd^2) %*% rep(1, ns)
    w <- w * (matrix(rep(s2, ns), ncol = ns) - wd * matrix(rep(s1,
        ns), ncol = ns))
    w <- w/(matrix(rep(s2, ns), ncol = ns) * matrix(rep(s0, ns),
        ncol = ns) - matrix(rep(s1, ns), ncol = ns)^2)
    estimate <- w %*% yw
    if (!(opt$display %in% "none"))
        lines(eval.points, estimate, lty = opt$lty)
    invisible(list(estimate = estimate, eval.points = eval.points,
        h = h, hv = hv, call = match.call()))
}
"sm.ts.pdf" <-
function (x, h = hnorm(x), lags, maxlag = 1, ask = TRUE)
{
    if (missing(lags))
        lags <- (1:maxlag)
    else maxlag <- max(lags)
    if (any(diff(lags)) < 0)
        stop("lags must be in increasing order")
    x.name <- deparse(substitute(x))
    x <- as.vector(x)
    n <- length(x)
    marginal <- sm.density(x, ylab = "Marginal density", xlab = x.name)
    if (ask)
        pause()
    for (m in lags) {
        x1 <- x[(m + 1):n]
        x0 <- x[1:(n - m)]
        biv <- sm.density(cbind(x0, x1), h = rep(h, 2), xlab = paste(x.name,
            "(t-", as.character(m), ")", sep = ""), ylab = paste(x.name,
            "(t)", sep = ""))
        biv$lag <- m
        title(paste("Density of lagged data of ", x.name, " (lag=",
            as.character(m), ")", sep = ""))
        if (ask & (m < maxlag))
            pause()
    }
    invisible(list(marginal = marginal, bivariate = biv))
}
"sm.weight" <-
function (x, eval.points, h, cross = FALSE, weights = rep(1, n), options)
{
    if (!exists(".sm.Options"))
        stop("cannot find .sm.Options")
    opt <- sm.options(options)
    replace.na(opt, hmult, 1)
    replace.na(opt, h.weights, rep(1, length(x)))
    replace.na(opt, poly.index, 1)
    poly.index <- opt$poly.index
    h.weights <- opt$h.weights
    hmult <- opt$hmult
    n <- length(x)
    ne <- length(eval.points)
    wd <- matrix(rep(eval.points, rep(n, ne)), ncol = n, byrow = TRUE)
    wd <- wd - matrix(rep(x, ne), ncol = n, byrow = TRUE)
    w <- matrix(rep(h.weights, ne), ncol = n, byrow = TRUE)
    w <- exp(-0.5 * (wd/(h * hmult * w))^2)
    w <- w * matrix(rep(weights, ne), ncol = n, byrow = TRUE)
    if (cross)
        diag(w) <- 0
    if (poly.index == 0) {
        den <- w %*% rep(1, n)
        w <- w/matrix(rep(den, n), ncol = n)
    }
    else if (poly.index == 1) {
        s0 <- w %*% rep(1, n)
        s1 <- (w * wd) %*% rep(1, n)
        s2 <- (w * wd^2) %*% rep(1, n)
        w <- w * (matrix(rep(s2, n), ncol = n) - wd * matrix(rep(s1,
            n), ncol = n))
        w <- w/(matrix(rep(s2, n), ncol = n) * matrix(rep(s0,
            n), ncol = n) - matrix(rep(s1, n), ncol = n)^2)
    }
}
"sm.weight2" <-
function (x, eval.points, h, cross = FALSE, weights = rep(1, nrow(x)),
    options = list())
{
    opt <- sm.options(options)
    n <- nrow(x)
    ne <- nrow(eval.points)
    replace.na(opt, h.weights, rep(1, n))
    h.weights <- opt$h.weights
    hmult <- opt$hmult
    wd1 <- matrix(rep(eval.points[, 1], rep(n, ne)), ncol = n,
        byrow = TRUE)
    wd1 <- wd1 - matrix(rep(x[, 1], ne), ncol = n, byrow = TRUE)
    w <- matrix(rep(h.weights, ne), ncol = n, byrow = TRUE)
    w <- exp(-0.5 * (wd1/(h[1] * hmult * w))^2)
    wd2 <- matrix(rep(eval.points[, 2], rep(n, ne)), ncol = n,
        byrow = TRUE)
    wd2 <- wd2 - matrix(rep(x[, 2], ne), ncol = n, byrow = TRUE)
    w <- w * exp(-0.5 * (wd2/(h[2] * hmult * matrix(rep(h.weights,
        ne), ncol = n, byrow = TRUE)))^2)
    w <- w * matrix(rep(weights, ne), ncol = n, byrow = TRUE)
    if (cross)
        diag(w) <- 0
    if (opt$poly.index == 0) {
        den <- w %*% rep(1, n)
        w <- w/matrix(rep(den, n), ncol = n)
    }
    else if (opt$poly.index == 1) {
        a11 <- w %*% rep(1, n)
        a12 <- (w * wd1) %*% rep(1, n)
        a13 <- (w * wd2) %*% rep(1, n)
        a22 <- (w * wd1^2) %*% rep(1, n)
        a23 <- (w * wd1 * wd2) %*% rep(1, n)
        a33 <- (w * wd2^2) %*% rep(1, n)
        d <- a22 * a33 - a23^2
        b1 <- 1/(a11 - ((a12 * a33 - a13 * a23) * a12 + (a13 *
            a22 - a12 * a23) * a13)/d)
        b2 <- (a13 * a23 - a12 * a33) * b1/d
        b3 <- (a12 * a23 - a13 * a22) * b1/d
        wt <- matrix(rep(b1, n), ncol = n)
        wt <- wt + matrix(rep(b2, n), ncol = n) * wd1
        wt <- wt + matrix(rep(b3, n), ncol = n) * wd2
        w <- wt * w
    }
    w
}
"sphdraw" <-
function (theta, phi)
{
    a1 <- 0
    a2 <- 30
    a3 <- 60
    a4 <- 90
    a5 <- 120
    a6 <- 150
    b1 <- (-90)
    b2 <- (-60)
    b3 <- (-30)
    b4 <- 0
    b5 <- 30
    b6 <- 60
    b7 <- 90
    latlines(b1, theta, phi)
    latlines(b2, theta, phi)
    latlines(b3, theta, phi)
    latlines.e(b4, theta, phi)
    latlines(b5, theta, phi)
    latlines(b6, theta, phi)
    latlines(b7, theta, phi)
    longlines.e(a1, theta, phi)
    longlines(a2, theta, phi)
    longlines(a3, theta, phi)
    longlines(a4, theta, phi)
    longlines(a5, theta, phi)
    longlines(a6, theta, phi)
    circle(1)
}
"sphimage" <-
function (latitude, longitude, kap, theta, phi, ngrid = 32)
{
    values <- seq(-1 + 1/ngrid, 1 - 1/ngrid, length = ngrid)
    xgrid <- rep(values, rep(ngrid, ngrid))
    ygrid <- rep(values, ngrid)
    dvec <- rep(0, ngrid^2)
    xlong <- longitude * pi/180
    xlat <- latitude * pi/180
    n <- length(longitude)
    radtheta <- theta * pi/180
    radphi <- phi * pi/180
    xgrid[xgrid^2 + ygrid^2 >= 1] <- NA
    ygrid[xgrid^2 + ygrid^2 >= 1] <- NA
    za <- -xgrid * sin(radtheta) - ygrid * cos(radtheta) * sin(radphi)
    zb <- cos(radphi) * cos(radtheta) * sqrt(1 - xgrid^2 - ygrid^2)
    z <- za + zb
    if ((theta == 90) | (theta == 270))
        x <- -ygrid * sin(radtheta) * sin(radphi) + cos(radphi) *
            sqrt(1 - ygrid^2 - z^2)
    else x <- (xgrid + z * sin(radtheta))/cos(radtheta)
    if (phi == 90)
        y <- sqrt(1 - x^2 - z^2)
    else if (phi == -90)
        y <- -sqrt(1 - x^2 - z^2)
    else y <- (ygrid + (x * sin(radtheta) + z * cos(radtheta)) *
        sin(radphi))/cos(radphi)
    xyzok <- (((x/sqrt(x^2 + z^2)) * (sqrt(1 - y^2)) * sin(radtheta) *
        cos(radphi)) + (y * sin(radphi)) - ((-z/sqrt(x^2 + z^2)) *
        (sqrt(1 - y^2)) * cos(radphi) * cos(radtheta)))
    other <- !is.na(xyzok) & xyzok < 0
    z[other] <- (za - zb)[other]
    x[other] <- ((xgrid + (z * sin(radtheta)))/cos(radtheta))[other]
    y[other] <- ((ygrid + ((x * sin(radtheta)) + (z * cos(radtheta))) *
        sin(radphi))/cos(radphi))[other]
    xj <- cos(xlong) * cos(xlat)
    yj <- sin(xlat)
    zj <- -sin(xlong) * cos(xlat)
    dvec <- exp(kap * cbind(x, y, z) %*% rbind(xj, yj, zj)) %*% rep(1/n, n)
    dvec[is.na(xgrid)] <- 0
    dvec <- dvec/max(dvec)
    fmat <<- matrix(dvec, ngrid, ngrid, byrow = TRUE)
    x <- seq(-1 + 1/ngrid, 1 - 1/ngrid, length = ngrid)
    y <- x
    image(x, y, fmat, add = TRUE)
    angle <- seq(0, pi/2, length = 50)
    xx <- cos(angle)
    yy <- sin(angle)
    polygon(c(xx, 0, 1, 1), c(yy, 1, 1, 0), col = 0, border = 0)
    angle <- seq(pi/2, pi, length = 50)
    xx <- cos(angle)
    yy <- sin(angle)
    polygon(c(xx, -1, -1, 0), c(yy, 0, 1, 1), col = 0, border = 0)
    angle <- seq(pi, 3 * pi/2, length = 50)
    xx <- cos(angle)
    yy <- sin(angle)
    polygon(c(xx, 0, -1, -1), c(yy, -1, -1, 0), col = 0, border = 0)
    angle <- seq(3 * pi/2, 2 * pi, length = 50)
    xx <- cos(angle)
    yy <- sin(angle)
    polygon(c(xx, 1, 1, 0), c(yy, 0, -1, -1), col = 0, border = 0)
    sphdraw(theta, phi)
}
"type" <- function (descr = "", x, digits = 4)
{
    cat(paste(descr, " "))
    cat(paste(round(x, digits = digits)))
    cat("\n")
}
"wmean" <- function (x, w)
    sum(x * w)/sum(w)
"wvar" <- function (x, w)
    sum(w * (x - wmean(x, w))^2)/(sum(w) - 1)
