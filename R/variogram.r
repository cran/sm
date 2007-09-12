"sm.variogram" <- function(x, y, h, ...) {

   data    <- sm.check.data(x = x, y = y, ...)
   x       <- data$x
   y       <- data$y
   n       <- data$nobs
   ndim    <- data$ndim
   opt     <- data$options
   rawdata <- list(x = x, y = y, nbins = opt$nbins, nobs = n, ndim = ndim)

   replace.na(opt, display, "means")
   replace.na(opt, ngrid,   100)
   replace.na(opt, band,    TRUE)
   replace.na(opt, test,    TRUE)
   replace.na(opt, col,     "black")
   replace.na(opt, se,      FALSE)
   replace.na(opt, df,      5)
   replace.na(opt, nbins,   100)
   if (opt$nbins == 0) opt$nbins <- 100
   if (ndim == 1)
      x <- matrix(x, ncol = 1)

   x1mat <- matrix(rep(x[, 1], n), ncol = n, byrow = TRUE)
   if (ndim == 2) {
      x2mat <- matrix(rep(x[, 2], n), ncol = n)
      hmat  <- sqrt((t(x1mat) - x1mat)^2 + (t(x2mat) - x2mat)^2)
      }
   else
      hmat <- abs(t(x1mat) - x1mat)
   dmat <- matrix(rep(y, n), ncol = n)
   dmat <- sqrt(abs(t(dmat) - dmat))
   imat <- matrix(rep(1:n, n), ncol = n, byrow = TRUE)
   ind  <- as.vector(t(imat) - imat)
   hall <- (as.vector(hmat))[ind > 0]
   dall <- (as.vector(dmat))[ind > 0]
   i1   <- (as.vector(imat))[ind > 0]
   i2   <- (as.vector(t(imat)))[ind > 0]

#---------Bin the differences-------------

   bins  <- binning(hall, dall, nbins = opt$nbins)
   hh    <- bins$x
   dd    <- bins$means
   wts   <- bins$x.freq
   brks  <- bins$breaks
   empse <- sqrt(bins$devs / (wts - 1)) / sqrt(wts)
   empse[wts == 1] <- 0
   igp   <- as.vector(cut(hall, brks, labels = FALSE))

   if(missing(h)) h <- h.select(hh, dd, weights = wts, ...)
   b <- h

#---Construct covariance matrix of binned data---

   vv    <- 0.1724
   cv    <- 0.03144
   Sigma <- table(c(igp, igp), c(i1, i2))
   Sigma <- Sigma %*% t(Sigma)
   Sigma <- cv * (Sigma - diag(2 * wts))
   Sigma <- Sigma / outer(wts, wts)
   Sigma <- diag(vv / wts) + Sigma

#--------Implement test---------

   if (opt$test) {
      W    <- sm.weight(hh, hh, b, weights = wts, options = opt)
      est  <- W %*% dd
      r0   <- sum(wts * (dd - mean(dall))^2)
      r1   <- sum(wts * (dd - est)^2)
      tobs <- (r0 - r1) / r1
      nb   <- length(hh)
      A    <- matrix(rep(wts / sum(wts), nb), ncol = nb, byrow = TRUE)
      A    <- t(diag(nb) - A) %*% diag(wts) %*% (diag(nb) - A)
      A    <- A - (1 + tobs) * t(diag(nb) - W) %*% diag(wts) %*% (diag(nb) - W)
      pval <- p.quad.moment(A, Sigma, 0, 0)
      if (opt$verbose > 0) 
         cat("Test of spatial independence: p = ",round(pval, 3), "\n")
      }

#-------------Plot the data---------------

   replace.na(opt, eval.points, seq(min(hall), max(hall), length = opt$ngrid))
   ev  <- opt$eval.points
   W   <- sm.weight(hh, ev, b, weights = wts, options = opt)
   est <- W %*% dd
      
   if (opt$band | opt$se) {
      sigmahat <- sqrt(var(y))
      nmeans   <- length(wts)
      V        <- matrix(rep(wts / sum(wts), length(ev)), ncol = nmeans, byrow = TRUE)
      se.band  <- sigmahat * sqrt(diag((W - V) %*% Sigma %*% t(W - V)))
      se       <- sigmahat * sqrt(diag(Sigma))
      }

   if (opt$display != "none") {
   	
   	  if (!opt$add) {
         if (opt$display == "means") {
         	 xx <- hh
         	 yy <- dd
            }
         else {
         	 xx <- hall
         	 yy <- dall
            }
      	 replace.na(opt, xlab, "Distance")
      	 replace.na(opt, ylab, "Square-root difference")
         replace.na(opt, xlim, range(xx))
         if (opt$se & opt$band) {
            r <- range(yy, dd - 2 * se, dd + 2 * se, 
                           mean(dall) + 2 * se.band, mean(dall) - 2 * se.band)
            replace.na(opt, ylim, r)
            }
         else if (opt$se)
            replace.na(opt, ylim, range(yy, dd - 2 * se, dd + 2 * se))
         else if (opt$band)
            replace.na(opt, ylim, range(yy, mean(dall) + 2 * se.band, mean(dall) - 2 * se.band))
         else
            replace.na(opt, ylim, range(yy))
         plot(xx, yy, xlab = opt$xlab, ylab = opt$ylab, xlim = opt$xlim, ylim = opt$ylim, type = "n")
   	     }

      if (opt$band)
         polygon(c(ev, rev(ev)), c(mean(dall) + 2 * se.band, rev(mean(dall) - 2 * se.band)),
                 border = FALSE, col = opt$col.band)

      if (opt$display == "means") {
         points(hh, dd, col = opt$col.points, pch = opt$pch)
         if (opt$se) segments(hh, dd - 2 * se, hh, dd + 2 * se, col = opt$col.points)
         }
      else
         points(h, dall, col = opt$col.points, pch = opt$pch)

      lines(ev, est, col = opt$col, lty = opt$lty)
      }

#-------------Return values---------------

   results <- list(distance = hall, sqrtdiff = dall, eval.points = ev, estimate = est,
                   distance.mean = hh, sqrtdiff.mean = dd, weights = wts, h = h, 
                   ibin = match(igp, sort(unique(igp))), ipair = cbind(i1, i2))
   if (opt$test) results$p       <- pval
   if (opt$se)   results$se      <- se
   if (opt$band) results$se.band <- se.band
   invisible(results)
   
   }
