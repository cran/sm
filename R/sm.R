# This is the S-plus `sm library' for nonparametric smoothing methods
# Copyright (C) 1997 A.W.Bowman & A.Azzalini
# Update: 20 Aug 1997

sm.ancova <- function(x, y, group, h,
			display = "lines",
			model = "none", band = T, test = T,
			h.alpha =  2 * diff(range(x)) / length(x),
			ngrid = 50, eval.points = NA,
			xlab, ylab, ...)

{

if (model=="none") {
	band <- F
	test <- F
	}
if (missing(xlab))  xlab <- deparse(substitute(x))
if (missing(ylab))  ylab <- deparse(substitute(y))

fact        <- factor(group)
fact.levels <- levels(fact)
nlevels     <- length(fact.levels)

if (any(is.na(eval.points))) {
  start.eval <- max(tapply(x, fact, min))
  stop.eval  <- min(tapply(x, fact, max))
  eval.points <- seq(start.eval, stop.eval, length = ngrid)
  }

if (!(display == "none")) {
  plot(x, y, type = "n", xlab = xlab, ylab = ylab, ...)
  text(x, y, as.character(fact))

  if (!band) {
  for (i in 1:nlevels) {
    sm.regression(x[fact == fact.levels[i]], y[fact == fact.levels[i]],
	h = h, ngrid = ngrid, add = T, lty = i)
      }
    }
  }

ord <- order(fact, x)
xx  <- x[ord]
yy  <- y[ord]
fac <- fact[ord]
fac.levels <- levels(fac)
n   <- table(fac)

B  <- diag(0, sum(n))
Sd <- diag(0, sum(n))
istart <- 1
for (i in 1:nlevels) {
  xi <- xx[fac == fac.levels[i]]
  B[istart:(istart + n[i] - 1),
    istart:(istart + n[i] - 1)]  <- sm.sigweight(xi)
  Sd[istart:(istart + n[i] - 1),
    istart:(istart + n[i] - 1)]  <- sm.weight(xi, xi, h)
  istart <- istart + n[i]
  }
Ss <- sm.weight(xx, xx, h)
B <- B / (sum(n) - 2*nlevels)
sigma <- sqrt((yy %*% B %*% yy)[1,1])

if (model == "equal") {
  Q <- Sd - Ss
  Q <- t(Q) %*% Q
  obs <- ((yy %*% Q %*% yy) / sigma^2)[1,1]
  covar <- diag(sum(n))
  }

if (model == "parallel") {
  D <- matrix(0, ncol = nlevels - 1, nrow = sum(n))
  istart <- n[1] + 1
  for (i in 2:nlevels) {
    	D[istart:(istart + n[i] - 1),i - 1] <- 1
    }
  Q <- diag(sum(n)) - sm.weight(xx, xx, h.alpha)
  Q <- solve(t(D) %*% t(Q) %*% Q %*% D) %*% t(D) %*% Q %*% Q
  alpha <- as.vector(Q %*% yy)
  covar <- rbind( cbind(Q %*% t(Q), Q), cbind(t(Q), diag(sum(n))))
  adjyy <- yy - D %*% alpha
  ghat  <- Ss %*% adjyy
  ghati <- Sd %*% yy
  obs   <- sum((as.vector(D %*% alpha) + ghat - ghati)^2) / sigma^2
  Q     <- cbind((diag(sum(n)) - Ss) %*% D, (Ss - Sd))
  Q     <- t(Q) %*% Q
  B     <- rbind(matrix(0, nrow = nlevels-1, ncol = sum(n)+nlevels-1),
                 cbind(matrix(0, nrow = sum(n), ncol = nlevels-1), B)
                 )
  }

p <- NULL
if (!(model == "none")) {
  p   <- p.quad.moment(Q - B * obs, covar)
  cat("Test of ", model, " lines:   h = ",
	signif(h), "   p-value = ", round(p, 4), "\n")
  }

if (band & !(display == "none")) {
  if (nlevels > 2) print("Band available only to compare two groups.")
  else {
    model1 <- sm.regression(x[fact == fact.levels[1]],
 	y[fact == fact.levels[1]], h = h, eval.points = eval.points,
 	display = "none", ngrid = ngrid, add = T, lty = 1)
    model2 <- sm.regression(x[fact == fact.levels[2]],
 	y[fact == fact.levels[2]], h = h, eval.points = eval.points,
 	display = "none", ngrid = ngrid, add = T, lty = 2)
    model.y <- (model1$estimate + model2$estimate) / 2
    if (model == "parallel"){
       model.y <- cbind(model.y - alpha/2, model.y + alpha/2)
       }
    se <- sqrt((model1$se/model1$sigma)^2 + (model2$se/model2$sigma)^2)
    se <- se * sigma
    upper <- model.y + se
    lower <- model.y - se
    if (model == "equal") {
      upper <- pmin(pmax(upper, par()$usr[3]), par()$usr[4])
      lower <- pmin(pmax(lower, par()$usr[3]), par()$usr[4])
      polygon(c(eval.points, rev(eval.points)), c(lower, rev(upper)),
		border = F, col = "cyan")
      }
    else if (model == "parallel") {
      upper[,1] <- pmin(pmax(upper[,1], par()$usr[3]), par()$usr[4])
      lower[,1] <- pmin(pmax(lower[,1], par()$usr[3]), par()$usr[4])
      upper[,2] <- pmin(pmax(upper[,2], par()$usr[3]), par()$usr[4])
      lower[,2] <- pmin(pmax(lower[,2], par()$usr[3]), par()$usr[4])
      polygon(c(eval.points, rev(eval.points)),
              c(lower[,1],   rev(upper[,1])), col=6)
      polygon(c(eval.points, rev(eval.points)),
              c(lower[,2],   rev(upper[,2])), col=5)
      }
    text(x, y, as.character(fact))
    lines(eval.points, model1$estimate, lty = 1)
    lines(eval.points, model2$estimate, lty = 2)
    }
  }

r <- list(p = p, model = model, sigma = sigma)
if (model == "parallel") r <- list(p = p, model = model, sigma = sigma,
					alphahat = alpha)
invisible(r)

}


sm.sigweight <- function(x) {

	n   <- length(x)
	xx  <- sort(x)
	xx1 <- diff(xx)
	xx2 <- diff(xx, lag = 2)

	a   <- xx1[-1]/xx2
	b   <- xx1[-(n-1)]/xx2
        a[xx2==0] <- 0.5
        b[xx2==0] <- 0.5
	c   <- sqrt(a^2 + b^2 + 1)

	D   <- cbind(rep(0,n-2), diag(-1/c), rep(0,n-2)) +
		cbind(diag(a/c), rep(0,n-2), rep(0,n-2)) +
		cbind(rep(0,n-2), rep(0,n-2), diag(b/c))
	D   <- rbind(rep(0,n), D, rep(0,n))

	t(D) %*% D

}

binning <- function(x, breaks, nbins)
{# converts x to a frequency table; switches to 1-d or 2-d case
   if(is.vector(x) | (isMatrix(x) & min(dim(x))==1))
      {# 1-d case
      x <- as.vector(x)
      if(missing(nbins)) nbins <- round(log(length(x),2)+1)
      if(missing(breaks))  {
         breaks<-seq(min(x),max(x),length=nbins+1)
         breaks[1] <- breaks[1]-10^(-5)
         }
      else nbins<-length(breaks)-1
      result<-binning1(x,breaks=breaks,nbins=nbins)
      }
   else # otherwise, assume 2-d case
      {
      if(!isMatrix(x)) stop("wrong parameter x for binning")
      if(missing(nbins)) nbins <- round(log(dim(x)[1],2)+1)
      if(missing(breaks)) {
         breaks<-cbind( seq(min(x[,1]),max(x[,1]),length=nbins+1),
                        seq(min(x[,2]),max(x[,2]),length=nbins+1))
         breaks[1,]<-breaks[1,]-rep(10^(-5),2)
         }
      else nbins<-length(breaks)/2-1
      result<-binning2(x,breaks=breaks,nbins=nbins)
      }
   result
}

binning1 <- function(x, breaks, nbins=round(log(length(x),2)+1))
{# 1-d binning: output is vector of midpoints of bins, and vector of frequencies
   f<-cut(x,breaks=breaks)  # include.lowest=T)
   if(any(is.na(f)))  stop("breaks do not span the range of x")
   freq<-tabulate(f,length(levels(f)))
   midpoints<-(breaks[-1]+breaks[-(nbins+1)])/2
   list(x=midpoints, freq=freq)
}


binning2<-function(x, breaks, nbins=round(log(length(x)/2,2)+1))
{# 2-d binning: output is (nbins x 2) matrix of midpoints on the two axes,
 #              and a (nbins x nbins) matrix of frequencies
   f1<-cut(x[,1],breaks=breaks[,1]) # include.lowest=T)
   f2<-cut(x[,2],breaks=breaks[,2]) # include.lowest=T)
   freq<-t(table(f1,f2))
   dimnames(freq)<-NULL
   midpoints<-(breaks[-1,]+breaks[-(nbins+1),])/2
   z1 <- midpoints[,1]
   z2 <- midpoints[,2]
   X  <- cbind(rep(z1,length(z2)), rep(z2,rep(length(z1),length(z2))))
   X.f<- as.vector(t(freq))
   id <-(X.f>0)
   X  <-X[id,]
   X.f<-X.f[id]
   list(x=X, x.freq=X.f, midpoints=midpoints, breaks=breaks, freq.table=freq)
}


britmap <- function() {

#..................Plot map of British Isles................................

   provide.data(britpts)
   jump <- c(0, sqrt(diff(britlat)^2 + diff(britlong)^2))
   flag <- rep(1, nrow(britpts))
   flag[jump>=0.6] <- NA
   lines(britpts * flag)

  }

#            Smooth density estimation in  1 or 2 dimensions

sm.density <- function(x, h, hmult = 1,
	h.weights = NA,
	model = "none", band = T, add = F, lty = 1,
	display , props = c(75,50,25),
	xlab = NA, ylab = NA, zlab = NA,
	xlim = NA, ylim = NA, yht = NA, ngrid = NA, eval.points = NA,
	panel = F, positive=F, delta, weights=rep(1,n), theta, phi, ...)
{
xlab0 <- deparse(substitute(x)); xlab0

if (length(dim(x))>0) {
     ndim <- dim(x)[2]
     if (ndim>3) ndim <- 3
     n<-dim(x)[1]
     }
  else
     { ndim <- 1; x <- as.vector(x); n<-length(x)}

if(missing(delta)) {if(ndim==1) delta<-min(x) else delta<-apply(x,2,min)}

if(missing(h)){
   if(positive) {
     if(is.vector(x)) h <-hnorm(log(x+delta))
                else h <-hnorm(log(x+outer(rep(1,n),delta)))
     }
   else h<-hnorm(x)}

if(missing(display)) {if (positive) display<-"slice" else display<-"persp"}

if (ndim==1) {
	if (is.na(xlab)) xlab <- xlab0
        est <- sm.density.1d(x, h, hmult, h.weights, kernel.index=2,
	          model, band, display, add, lty, panel,
		  xlab, xlim, ylim, yht, ngrid, eval.points,
                  positive, delta, weights, ...)
	}

if (ndim==2) {
	if (is.na(xlab)) {
		if(!is.null(attributes(x)$dimnames))
		   xlab <- attributes(x)$dimnames[[2]][1]
		else xlab <- paste(xlab0,"[1]")
		}
	if (is.na(ylab)) {
		if(!is.null(attributes(x)$dimnames))
		   ylab <- attributes(x)$dimnames[[2]][2]
		else ylab <- paste(xlab0,"[2]")
		}
        if(missing(delta)) delta<-c(min(x[,1]),min(x[,2]))
	est <- sm.density.2d(x, h, hmult, h.weights, display, props,
		add, panel,
		xlim, ylim, xlab, ylab, zlab, ngrid, eval.points,
                positive=positive, delta=delta, lty, weights, ...)
	}

if (ndim==3) {
	if (is.na(xlab)) {
		if(!is.null(attributes(x)$dimnames))
		   xlab <- attributes(x)$dimnames[[2]][1]
		else if (!is.null(attributes(x)$names))
		   xlab <- attributes(x)$names[1]
		else xlab <- paste(xlab0,"[1]")
		}
	if (is.na(ylab)) {
		if(!is.null(attributes(x)$dimnames))
		   ylab <- attributes(x)$dimnames[[2]][2]
		else if (!is.null(attributes(x)$names))
		   ylab <- attributes(x)$names[2]
		else ylab <- paste(xlab0,"[2]")
		}
	if (is.na(zlab)) {
		if(!is.null(attributes(x)$dimnames))
		   zlab <- attributes(x)$dimnames[[2]][3]
		else if (!is.null(attributes(x)$names))
		   zlab <- attributes(x)$names[3]
		else zlab <- paste(xlab0,"[3]")
		}
	if (missing(theta)) theta <- pi/4
	if (missing(phi))   phi   <- pi/4
	est <- sm.density.3d(x, h, hmult, contour=props[1], ngrid=ngrid,
		xlab=xlab, ylab=ylab, zlab=zlab, theta=theta, phi=phi)
	}

invisible(est)

}
#             Smooth density estimation in 1 dimension

sm.density.1d <- function(x, h=hnorm(x), hmult = 1,
		   h.weights = rep(1,length(x)), kernel.index = 2,
		   model = "none", band = T,
		   display = "estimate", add = F, lty = 1, panel = F,
		   xlab = NA, xlim , ylim, yht, ngrid , eval.points = NA,
                   positive=F, delta, weights, ...)
{


#................Check and set function arguments......................

absent<-function(x) missing(x) | any(is.na(x))

if (absent(h.weights)) h.weights <- rep(1,length(x))
      else  band <- panel <- F
if (model=="none") band <- F
if (add | display=="none") panel <- F
if (absent(xlim)){
      if(positive) xlim<-c(0,max(x))
      else xlim <- c(min(x) - diff(range(x))/4, max(x) + diff(range(x))/4)
    }
if(absent(yht)){
   if (positive)
      yht<- max(0.4/(quantile(x,0.75)-quantile(x,0.50)),
                0.4/(quantile(x,0.50)-quantile(x,0.25)))
   else
      yht<-0.6/sqrt(var(x))
   }

if (absent(ylim)) ylim <- c(0, yht)

if (absent(ngrid)) ngrid <- 100

if (!add & !(display=="none"))
	plot(xlim, ylim, type="n", xlab = xlab, ylab = "Density",...)


#............................Plot density.................................

if (!(display=="none")) est <- smplot.density(x,h,hmult,h.weights,kernel.index,
				   band,add,lty,xlim,yht,ngrid,display,
                                   positive, delta, weights, ...)


#............................Panel control................................

if (panel) {
	items <-      c("Bandwidth:",
			"  - normal optimal",
			"  - plug-in",
			"  - xval",
			"  - increase",
			"  - decrease",
			"  - movie up",
			"  - movie down",
                        "Add Normal band",
			"Exit")
	hsj.flag <- F
	hcv.flag <- F

	ind <- menu(items, graphics=T, title="Density estimation")

	while (items[ind]!="Exit") {
		if (items[ind]=="  - normal optimal") {
			h <- hnorm(x)
			hmult <- 1
			}
		if (items[ind]=="  - plug-in") {
			if (!hsj.flag) {
				h.sj     <- hsj(x)
				hsj.flag <- T
				}
			h <- h.sj
			hmult <- 1
			}
		if (items[ind]=="  - xval") {
			if (!hcv.flag) {
				h.cv     <- hcv(x)
				hcv.flag <- T
				}
			h <- h.cv
			hmult <- 1
			}
		else if (items[ind]=="  - increase") {
			hmult <- hmult * 1.1
			}
		else if (items[ind]=="  - decrease") {
			hmult <- hmult / 1.1
			}
		else if (items[ind]=="  - movie up") {
			for (i in 1:6) {
			   hmult <- hmult * 1.1
			   est <- smplot.density(x,h,hmult,h.weights,kernel.index,
				band,add,lty,xlim,yht,ngrid,display,...)
			   }
			hmult <- hmult * 1.1
			}
		else if (items[ind]=="  - movie down") {
			for (i in 1:6) {
			   hmult <- hmult / 1.1
			   est <- smplot.density(x,h,hmult,h.weights,kernel.index,
				band,add,lty,xlim,yht,ngrid,display,...)
			   }
			hmult <- hmult / 1.1
			}
		else if (items[ind]=="Add Normal band" |
				items[ind]=="Remove Normal band") {
			band <- !band
			if (items[ind]=="Add Normal band") {
				items[ind] <- "Remove Normal band"  }
			else (items[ind] <- "Add Normal band")
			}

		est <- smplot.density(x,h,hmult,h.weights,kernel.index,
			band,add,lty,xlim,yht,ngrid,display,...)
		cat("h = ", signif(h*hmult,7), "\n")

		ind <- menu(items, graphics=T, title="Density estimation")
		}
	}


#.........................Assign returned values............................


if (!absent(eval.points)) {
   if (positive) est <- sm.density.positive.1d(x,h,hmult,h.weights,
                   xlim, ngrid, eval.points, delta, weights)
   else est <- sm.density.eval.1d(x, h, hmult, h.weights, kernel.index,
			xlim, ngrid, eval.points, weights)
   }
else if(display=="none") {
   if (positive) est <- sm.density.eval.1d(x, h, hmult, h.weights,
   			kernel.index, xlim, ngrid, weights=weights)
   else est <- sm.density.eval.1d(x, h, hmult, h.weights, kernel.index,
			xlim, ngrid, eval.points=NA, weights)
   }

n <- length(x)
if (all(weights==rep(1,n)) & all(h.weights==rep(1,n))
			& kernel.index == 2 & positive==F) {
  se <- sqrt(dnorm(0, sd = sqrt(2)) / (4 * n * h))
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

#---------------------------------------------------------------smplot.density

smplot.density <- function(x, h, hmult, h.weights, kernel.index,
                   band, add, lty, xlim, yht, ngrid, display,
                   positive=F, delta, weights=rep(1,length(x)), ...)
{
   if(positive)
     {est <- sm.density.positive.1d(x,h,hmult,h.weights,
                   xlim=xlim, ngrid, eval.points=NA, delta, weights, ...)
      }
   else{
     est <- sm.density.eval.1d(x, h, hmult, h.weights, kernel.index,
                          xlim, ngrid, eval.points=NA, weights=weights, ...)
     if (band) normdens.band(x, h, hmult, xlim, yht, ngrid)
	 else  if (!add) polygon(c(par()$usr[1:2],par()$usr[2:1]),
           rep(c(par()$usr[3], par()$usr[4]*0.999),c(2,2)), col=0, border=F)
     }
   box()
   lines(est$eval.points, est$estimate, lty=lty, ...)
   n <- length(x)
   if (display == "se" & all(weights==rep(1,n)) & all(h.weights==rep(1,n))
			& kernel.index == 2) {
      se <- sqrt(dnorm(0, sd = sqrt(2)) / (4 * n * h))
      upper <- sqrt(est$estimate) + 2 * se
      lower <- pmax(sqrt(est$estimate) - 2 * se, 0)
      upper <- upper^2
      lower <- lower^2
      lines(est$eval.points, upper, lty = 3)
      lines(est$eval.points, lower, lty = 3)
      }
#  text(xlim[1], yht*0.95, paste("h = ", signif(h*hmult,3)),adj=0)

   invisible(est)
}


#--------------------------------------------------------- sm.density.eval.1d

sm.density.eval.1d <- function(x, h,
	  hmult=1, h.weights=rep(1,length(x)), kernel.index=2,
	  xlim=c(min(x) - diff(range(x))/4, max(x) + diff(range(x))/4),
	  ngrid=100, eval.points = NA, weights=rep(1,n), ...)
{
	if(any(is.na(eval.points)))
		xnew <- seq(xlim[1], xlim[2], length=ngrid)
	else xnew <- eval.points
	n     <- length(x)
	neval <- length(xnew)
	W     <- matrix(rep(xnew, rep(n, neval)), ncol = n, byrow = T)
	W     <- W - matrix(rep(x, neval), ncol = n, byrow = T)
	W1    <- matrix(rep(h.weights, neval), ncol = n, byrow = T)
	W     <- exp(-.5 * (W/(hmult*h*W1))^2) / W1
	est   <- W  %*% weights/(sum(weights)*sqrt(2*pi)*hmult*h)
	invisible(list(eval.points = xnew, estimate = as.vector(est),
			h=h*hmult, h.weights=h.weights, weights=weights))
}


sm.density.positive.1d<-function(x,h=hnorm(log(x+delta)),hmult=1,
   h.weights=rep(1,length(x)),
   ngrid=100, xlim=range(x), eval.points=NA, delta=min(x), ...)
{ # estimate density of positive variable x by transform to log(x+delta)
  if(min(x)<=0)   cat("Warning: some data are not positive\n")
  if(min(xlim)<0) cat("Warning: xlim<0 with positive=T \n")
  if(missing(eval.points) | any(is.na(eval.points))){
    a<-log(xlim+1/ngrid)
    eval.points<-exp(seq(min(a),max(a),length=ngrid))-1/ngrid
    }
  f<- sm.density.eval.1d(log(x+delta),h=h,hmult=1,h.weights,
         eval.points=log(eval.points+delta))
  est <- f$estimate/(eval.points+delta)
  est[is.na(est)]<-0
  list(eval.points=eval.points, estimate=as.vector(est),h=h)
}



#--------------------------------------------------------------normdens.band


normdens.band <- function(x, h, hmult, xlim, yht, ngrid) {
  x.points <- seq(xlim[1],xlim[2],length=ngrid)
  n        <- length(x)
  xbar     <- mean(x)
  sx       <- sqrt(var(x))
  hm       <- h*hmult
  dmean    <- dnorm(x.points,xbar,sqrt(sx^2+hm^2))
  dvar     <- (dnorm(0,0,sqrt(2*hm^2))*dnorm(x.points,xbar,sqrt(sx^2+0.5*hm^2))
                  -(dmean)^2)/n
  upper    <- pmin(dmean+2*sqrt(dvar), par()$usr[4])
  lower    <- pmax(0,dmean-2*sqrt(dvar))

#     The following line blanks out the current plotting area.
#     This makes animation a little smoother.
  polygon(c(par()$usr[1:2],par()$usr[2:1]),
           rep(c(par()$usr[3], par()$usr[4]*0.999),c(2,2)), col=0, border=F)
  polygon(c(x.points,rev(x.points)),c(upper,rev(lower)),col="cyan",border=F)
  }

#----------------------------------------------------------------------hnorm

hnorm <- function(x) {

   if(isMatrix(x)) {
      ndim <- ncol(x)
      n  <- nrow(x)
      sd <- sqrt(apply(x,2,var))
      if (ndim==2) hh <- sd*(1/n)^(1/6)
      if (ndim==3) hh <- sd * (4/(5 * n))^(1/7)
      hh
      }
   else {
      ndim <- 1
      n <- length(x)
      sd <- sqrt(var(x))
      sd*(4/(3*n))^(1/5)
      }

  }

#-----------------------Sheather-Jones bandwidth-----------------------------

sj <- function(x, h) {

  phi6 <- function(x) (x^6 - 15*x^4 + 45*x^2 - 15)*dnorm(x)
  phi4 <- function(x) (x^4 - 6*x^2 + 3)*dnorm(x)

  n <- length(x)
  lambda <- quantile(x,0.75)-quantile(x,0.25)
  a <- 0.920*lambda*n^(-1/7)
  b <- 0.912*lambda*n^(-1/9)
  W <- matrix(rep(x, rep(n, n)), ncol = n, byrow = T)
  W <- W - matrix(rep(x, n),     ncol = n, byrow = T)
  W1 <- matrix(phi6(W/b),ncol=n)
  tdb <- as.numeric(rep(1,n) %*% W1 %*% rep(1,n))
  tdb <- -tdb/(n*(n-1)*b^7)
  W1 <- matrix(phi4(W/a),ncol=n)
  sda <- as.numeric(rep(1,n) %*% W1 %*% rep(1,n))
  sda <- sda/(n*(n-1)*a^5)
  alpha2 <- 1.357 * (abs(sda/tdb))^(1/7) * h^(5/7)
  W1 <- matrix(phi4(W/alpha2),ncol=n)
  sdalpha2 <- as.numeric(rep(1,n) %*% W1 %*% rep(1,n))
  sdalpha2 <- sdalpha2/(n*(n-1)*alpha2^5)

  result <- (dnorm(0,sd=sqrt(2)) / (n * abs(sdalpha2)))^0.2 - h
  attributes(result)$names <- NULL
  as.double(result)

  }

hsj <- function(x) {

#	Sheather-Jones choice of bandwidth for univariate
#	density estimation.

#	Note that, for a more accurate value, the process could
#	be repeated from the current value, using hstep = 0.99 or 1.01.

  h0 <- hnorm(x)
  v0 <- sj(x, h0)
  if (v0 > 0) hstep <- 1.1 else hstep <- 0.9
  h1 <- h0 * hstep
  v1 <- sj(x, h1)
  while (v1*v0 > 0) {
    h0 <- h1
    v0 <- v1
    h1 <- h0 * hstep
    v1 <- sj(x, h1)
    }
  h0 + (h1 - h0) * abs(v0) / (abs(v0) + abs(v1))
  }


#	Cross-validation bandwidth

#	Cross-validation function for a 1-d density estimate
#	Note that in the 2-d case this function uses a single h which is
#	then scaled by the standard deviations of each variable.

cv <- function(x, h, h.weights=NA) {

   if (!isMatrix(x)) {

   	n     <- length(x)
	if (missing(h.weights) | (any(is.na(h.weights))))
		h.weights <- rep(1,n)
   	hcvff <- sum(dnorm(0,mean=0,sd=sqrt(2)*h*h.weights))/(n*(n-1))

	W     <- matrix(rep(x, rep(n, n)), ncol = n, byrow = T)
	W     <- W - matrix(rep(x, n),     ncol = n, byrow = T)
	W1    <- matrix(rep(h.weights^2, n), ncol = n, byrow = T)

	W2    <- exp(-.5 * (W/(h*sqrt(W1+t(W1))))^2) /
			(sqrt(2*pi)*h*sqrt(W1+t(W1)))
	hcvff <- hcvff + (sum(W2) - sum(diag(W2)))*(n-2)/(n*(n-1)^2)

	W2    <- exp(-.5 * (W/(h*sqrt(W1)))^2) / (sqrt(2*pi)*h*sqrt(W1))
	hcvff <- hcvff - (sum(W2) - sum(diag(W2)))*2/(n*(n-1))

   	}

   if (isMatrix(x)) {

	x1  <- x[,1]
	x2  <- x[,2]
	h1  <- h * sqrt(var(x1))
	h2  <- h * sqrt(var(x2))
   	n   <- length(x1)
	if (missing(h.weights) | (any(is.na(h.weights))))
		h.weights <- rep(1,n)

   	hcvff <- sum(dnorm(0,mean=0,sd=sqrt(2)*h1*h.weights)*
   			dnorm(0,mean=0,sd=sqrt(2)*h2*h.weights))/(n*(n-1))

	W     <- matrix(rep(x1, rep(n, n)),   ncol = n, byrow = T)
	W     <- W - matrix(rep(x1, n),       ncol = n, byrow = T)
	W1    <- matrix(rep(h.weights^2, n),  ncol = n, byrow = T)
	W2    <- exp(-.5 * (W/(h1*sqrt(W1+t(W1))))^2) /
			(sqrt(2*pi)*h1*sqrt(W1+t(W1)))
	W     <- matrix(rep(x2, rep(n, n)), ncol = n, byrow = T)
	W     <- W - matrix(rep(x2, n),     ncol = n, byrow = T)
	W2    <- W2 * exp(-.5 * (W/(h2*sqrt(W1+t(W1))))^2) /
			(sqrt(2*pi)*h2*sqrt(W1+t(W1)))
	hcvff <- hcvff + (sum(W2) - sum(diag(W2)))*(n-2)/(n*(n-1)^2)

	W2    <- exp(-.5 * (W/(h2*sqrt(W1)))^2) / (sqrt(2*pi)*h2*sqrt(W1))
	W     <- matrix(rep(x1, rep(n, n)), ncol = n, byrow = T)
	W     <- W - matrix(rep(x1, n),     ncol = n, byrow = T)
	W2    <- W2 * exp(-.5 * (W/(h1*sqrt(W1)))^2) / (sqrt(2*pi)*h1*sqrt(W1))
	hcvff <- hcvff - (sum(W2) - sum(diag(W2)))*2/(n*(n-1))

   	}

	hcvff
}


hcv <- function(x, y=NA, h.weights = NA, ngrid = 8, hstart = NA, hend = NA,
			display = "none", add = F, ...)
  {

#	Cross-validatory choice of bandwidth.
#	Note that in the 2-d case the function cv scales the single
#	value of h supplied.  hcv therefore scales the answer by the
#	sd's of each variable.

  if (length(dim(x))>0) {
    ndim <- 2
    n <- length(x[,1])
    }
  else {
    ndim <- 1
    n <- length(x)
    }

  if (is.na(hstart)) {
    if (ndim==1) hstart <- hnorm(x) / 10
    else if (any(is.na(y))) hstart <- hnorm(x[,1]/sqrt(var(x[,1]))) / 10
    	 else               hstart <- hnorm(x[,1]/sqrt(var(x[,1]))) / 4
    }
  if (is.na(hend)) {
    if (ndim == 1) hend <- hnorm(x) * 2
    else           hend <- hnorm(x[,1]/sqrt(var(x[,1]))) * 2
    }

  if (missing(h.weights) | (any(is.na(h.weights))))
	h.weights <- rep(1,n)

  cvgrid <- vector("numeric", length = ngrid)
  hgrid  <- log(hstart) + (log(hend) - log(hstart)) * (0:(ngrid-1))/(ngrid-1)
  if (any(is.na(y))) {
     for (i in 1:ngrid) cvgrid[i] <- cv(x, exp(hgrid[i]), h.weights)
     }
  else {
     if (ndim==1) for (i in 1:ngrid) {
        cvgrid[i] <- sum((y - sm.weight(x, x, exp(hgrid[i]),
        		h.weights=h.weights, cross=T) %*% y)^2)
        }
     if (ndim==2) for (i in 1:ngrid) {
        cvgrid[i] <- sum((y - sm.weight2(x, x,
           exp(hgrid[i]*c(sqrt(var(x[,1])), sqrt(var(x[,2])))),
        		h.weights=h.weights, cross=T) %*% y)^2)
        }
     }

  if (any(is.na(cvgrid))) {
    cat("\n")
    cat("hcv: some computations failed.","\n")
    cat("Try readjusting hstart and hend.", "\n")
    cat("hstart: ", hstart, "\n")
    cat("hend  : ", hend,   "\n")
    cat("\n")
    print(cbind(h=exp(hgrid), cv=cvgrid))
    stop()
    }

  ind    <- (1:ngrid)[cvgrid == min(cvgrid)]

  if (!(display=="none")) {
    if (!add) {
      if (display=="log") plot(hgrid,cvgrid, type="l",
      				xlab = "Log h", ylab = "CV", ...)
      else plot(exp(hgrid),cvgrid, type="l", xlab = "h", ylab = "CV", ...)
      }
    else {
      if (display=="log") lines(hgrid, cvgrid, ...)
      else lines(exp(hgrid), cvgrid, ...)
      }
    }

  if (ind == 1 | ind == ngrid) {
    cat("\n")
    cat("hcv: boundary of search area reached.","\n")
    cat("Try readjusting hstart and hend.", "\n")
    cat("hstart: ", hstart, "\n")
    cat("hend  : ", hend,   "\n")
    cat("\n")
    print(cbind(h=exp(hgrid), cv=cvgrid))
    stop()
    }

  v0 <- cvgrid[ind-1]
  v1 <- cvgrid[ind]
  v2 <- cvgrid[ind+1]
  l0 <- hgrid[ind-1]
  l1 <- hgrid[ind]
  l2 <- hgrid[ind+1]
  aa  <- (v1-v0-(l1-l0)*(v1-v2)/(l1-l2))/
  		(l1^2-l0^2-(l1^2-l2^2)*(l1-l0)/(l1-l2))
  bb  <- (v1-v2-aa*(l1^2-l2^2))/(l1-l2)
  cc  <- v0 - aa*l0^2 - bb*l0

  h <- exp(-bb/(2*aa))
  if (ndim == 1) result <- h
  else           result <- c(h*sqrt(var(x[,1])), h*sqrt(var(x[,2])))
  result

  }

nise <- function(y, hmult = 1) {

#	ISE between a density estimate constructed from standardised
#	data and a standard Normal distribution.

  y <- (y-mean(y))/sqrt(var(y))
  h <- hnorm(y) * hmult
  result <- dnorm(0,sd=sqrt(2+2*h^2))
  result <- result - 2*sm.density(y, h=sqrt(1+2*h^2), eval.points=0,
  			display="none")$estimate
  result <- result + mean(sm.density(y, h=sqrt(2)*h, eval.points=y,
  			display="none")$estimate)
  result
  }


nmise <- function(sd, n, h) {

#	MISE of a density estimate constructed from Normal data with
#	std.dev. sd

	dnorm(0,sd=sqrt(2)*h)/n +
		(1-1/n)*dnorm(0,sd=sqrt(2*(sd^2+h^2))) -
		2*dnorm(0,sd=sqrt(2*sd^2+h^2)) +
		dnorm(0,sd=sqrt(2)*sd)
	}
#	        Smooth density estimation in 2 dimensions

sm.density.2d <- function(X, h = hnorm(X), hmult = 1, h.weights = NA,
		display="persp",
		props=c(25,50,75), add = F, panel = F,
		xlim = range(X[,1]), ylim = range(X[,2]),
		xlab = NA, ylab = NA, zlab = "Density",
                ngrid = 50, eval.points = NA, positive=F, delta, lty=1,
                weights=rep(1,length(x)), ...)
{
#................Check and set function arguments......................

x <- X[,1]
y <- X[,2]

if (length(xlim)==1) xlim  <- range(X[,1])
if (length(ylim)==1) ylim  <- range(X[,2])
if (is.na(ngrid))    ngrid <- 50
if (display=="none") panel <- F
if (any(is.na(h.weights))) h.weights <- rep(1, length(x))


#............................Plot density.................................

if (display=="persp")
	est <- sm.persplot(x,y,h,hmult,h.weights,xlim,ylim,xlab,ylab,zlab,
                        ngrid,positive,delta,weights,...)
else if (display=="image")
	est <- sm.imageplot(x,y,h,hmult,h.weights,xlim,ylim,xlab,ylab,
                        ngrid,positive,delta,weights,... )
else if (display=="slice")
	est <- sm.sliceplot(x,y,h,hmult,h.weights,props,add,xlim,ylim,
                        xlab,ylab,ngrid, positive,delta,lty,weights,...)


#............................Panel control................................

if (panel) {
	items <-      c("Bandwidth:",
			"  - Normal optimal",
			"  - increase",
			"  - decrease",
			"Exit")

	ind <- menu(items, graphics=T, title="2-d density estimation")

	while (items[ind]!="Exit") {
		if (items[ind]=="  - Normal optimal") {
			hmult <- 1
			}
		else if (items[ind]=="  - increase") {
			hmult <- hmult * 1.1
			}
		else if (items[ind]=="  - decrease") {
			hmult <- hmult / 1.1
			}

	if (display=="persp")
		est <- sm.persplot(x,y,h,hmult,h.weights,xlim,ylim,xlab,ylab,
                                zlab,ngrid,positive,delta,weights,...)
	else if (display=="image")
		est <- sm.imageplot(x,y,h,hmult,h.weights,xlim,ylim,xlab,ylab,
                          ngrid,positive,delta,weights,...)
	else if (display=="slice") {
                est <- sm.sliceplot(x,y,h,hmult,h.weights,props,add,
                          xlim,ylim,xlab,ylab,ngrid,positive,delta,
                          lty,weights,...)
		}

		ind <- menu(items, graphics=T, title="2-d density estimation")
		}
	}


#.........................Assign returned values............................

if (isMatrix(eval.points))
   est <- sm.density.eval.2d(x,y,h,hmult,h.weights,xnew=eval.points[,1],
                ynew= eval.points[,2], eval.type="points",weights=weights)
else if(display=="none")
   est <- sm.density.eval.2d(x,y,h,hmult,h.weights,xlim,ylim,
                eval.type="grid", ngrid=ngrid, weights=weights)

n <- length(x)
if (all(weights==rep(1,n)) & all(h.weights==rep(1,n))) {
  se <- dnorm(0, sd = sqrt(2)) / sqrt(4 * n * h^2)
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

#-----------------------------Perspective plot-----------------------------


sm.persplot <- function(x, y, h=hnorm(cbind(x,y)),
		hmult = 1, h.weights = NA, xlim = range(x), ylim = range(y),
		xlab = "x1", ylab = "x2", zlab = "Density", ngrid=50,
                positive,delta,weights, ...)
{

#     Perspective of a nonparametric density estimate
#     for bivariate data.

     if (any(is.na(h.weights))) h.weights <- rep(1, length(x))
     if (is.na(zlab)) zlab <- "Density"

     xgrid <- seq(xlim[1],xlim[2],length=ngrid)
     ygrid <- seq(ylim[1],ylim[2],length=ngrid)
     if(!positive)
        dgrid <- sm.density.eval.2d(x,y,h,hmult,h.weights,xlim,ylim,xgrid,
                        ygrid,eval.type="grid",ngrid,weights)$estimate
     else{
        f<-sm.density.positive.grid(cbind(x,y),h,delta,ngrid,NA,xlim,ylim)
        xgrid<-f$eval.points[,1]
        ygrid<-f$eval.points[,2]
        dgrid<-f$estimate
        }
#     persp(xgrid, ygrid, dgrid, xlab=xlab, ylab=ylab, zlab=zlab, ...)
     persp(xgrid, ygrid, dgrid, theta = -30, phi = 40, d = 4)

     invisible(list(eval.points = cbind(xgrid, ygrid), estimate = dgrid,
                h=h*hmult, h.weights=h.weights, weights=weights))
}


#--------------------------------Image plot--------------------------------


sm.imageplot <-function(x, y, h,
		hmult = 1, h.weights = NA, xlim = range(x), ylim = range(y),
                xlab = "x1", ylab = "x2", ngrid=50, positive, delta, weights,
                ...)
{

#     Greyscale plot of a nonparametric density estimate
#     for bivariate data.

     if (any(is.na(h.weights))) h.weights <- rep(1, length(x))

     xgrid <- seq(xlim[1],xlim[2],length=ngrid)
     ygrid <- seq(ylim[1],ylim[2],length=ngrid)
     if(!positive)
        dgrid <- sm.density.eval.2d(x,y,h,hmult,h.weights,xlim,ylim,xgrid,
                        ygrid,eval.type="grid",ngrid,weights)$estimate
     else{
        f<-sm.density.positive.2d(cbind(x,y),h,delta,ngrid,NA,
                xlim,ylim,eval.type="grid")
        xgrid<-f$eval.points[,1]
        ygrid<-f$eval.points[,2]
        dgrid<-f$estimate
        }
     image(xgrid, ygrid, dgrid, xlab=xlab, ylab=ylab, ...)

     invisible(list(eval.points = cbind(xgrid, ygrid), estimate = dgrid,
                h=h*hmult, h.weights=h.weights, weights=weights))
}



#--------------------------------Slice plot--------------------------------

sm.sliceplot <- function(x, y, h,
		hmult = 1, h.weights = NA, props=c(25,50,75),  add = F,
		xlim = range(x), ylim = range(y),
                xlab = "x1", ylab = "x2", ngrid=50, positive,delta,lty,
                weights,...)
{

#     Percentile contours of a nonparametric density estimate
#     for bivariate data.

        if (any(is.na(h.weights))) h.weights <- rep(1, length(x))
        if (!add) plot(x,y,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,...)
        if(positive){
           f <- sm.density.positive.2d(cbind(x,y),h,delta,ngrid,
                   xlim=xlim,ylim=ylim,eval.type="grid")
           dgrid <- f$estimate
           xgrid <- f$eval.points[,1]
           ygrid <- f$eval.points[,2]
           }
        else{
           f <- sm.density.eval.2d(x,y,h,hmult,h.weights,xlim,ylim,ngrid=ngrid,
                                    eval.type="grid",weights=weights)
           dgrid <- f$estimate
           xgrid <- f$eval.points[,1]
           ygrid <- f$eval.points[,2]
           }
#...........Evaluate the density estimate at the observations..............

        if(positive)
           dobs <- sm.density.positive.2d(cbind(x,y),h,delta,ngrid,
                        cbind(x,y),xlim,ylim, eval.type="points")$estimate
        else
           dobs<-sm.density.eval.2d(x,y,h,hmult,h.weights,xnew=x,ynew=y,
                        weights=weights)$estimate

#..........Scale heights to values of props and draw each contour..........

        hts   <- quantile(rep(dobs,weights),prob=(100-props)/100)
        for(i in 1:length(props)) {
	   scale <- props[i]/hts[i]
           contour(xgrid, ygrid, dgrid*scale, level=hts[i]*scale, add=T,
           		lty=lty, ...)
	   }
	invisible(list(eval.points = cbind(xgrid, ygrid), estimate = dgrid,
                        h=h*hmult, h.weights=h.weights, weights=weights))

}


#---------------- merge   sm.density.grid & sm.density.pts

sm.density.eval.2d <- function(x, y, h, hmult=1, h.weights = NA, xlim, ylim,
                xnew, ynew, eval.type="points", ngrid=50, weights=rep(1,n))
{
        if(missing(xlim)) xlim<-range(x)
        if(missing(ylim)) ylim<-range(y)
        if(missing(xnew)) xnew<-seq(xlim[1],xlim[2],length=ngrid)
        if(missing(ynew)) ynew<-seq(ylim[1],ylim[2],length=ngrid)
	n     <- length(x)
        if (any(is.na(h.weights))) h.weights <- rep(1, length(x))
	nnew  <- length(xnew)
	W1    <- matrix(rep(xnew, rep(n, nnew)), ncol = n, byrow = T)
	W1    <- W1 - matrix(rep(x   , nnew),    ncol = n, byrow = T)
	W2    <- matrix(rep(h.weights, nnew),    ncol = n, byrow = T)
	Wx    <- exp(-.5 * (W1/(hmult*h[1]*W2))^2) / W2
	W1    <- matrix(rep(ynew, rep(n, nnew)), ncol = n, byrow = T)
        W1    <- W1 - matrix(rep(y,    nnew),    ncol = n, byrow = T)
        Wy    <- exp(-.5 * (W1/(hmult*h[2]*W2))^2) / W2
	if(eval.type=="points")
           est <- ((Wx * Wy) %*% weights)/(sum(weights)*2*pi*h[1]*h[2]*hmult^2)
        else # (should be eval.type=="grid" then)
           est <- (Wx %*%(weights*t(Wy)))/(sum(weights)*2*pi*h[1]*h[2]*hmult^2)
	invisible(list(eval.points = cbind(xnew, ynew), estimate = est,
                        h=h*hmult, h.weights=h.weights, weights=weights))
}

#------------ case with (x,y) both positive rv's ------------------------

sm.density.positive.2d<-function(X,
       h=c(hnorm(log(X[,1]+delta[1])),hnorm(log(X[,2]+delta[2]))),
       delta=c(min(X[,1]),min(X[,2])), ngrid=50, eval.points=NA,
       xlim=range(X[,1]), ylim=range(X[,2]), eval.type="points",...)
{# estimate density of positive variables (x1,x2) by transform to
 # (log(x1+delta1), log(x2+delta2))

  if(min(X)<=0) cat("Warning: some data are not positive\n")
  if(dim(X)[2]!=2) cat("Data matrix does not have two columns\n")
  x1<-X[,1]
  x2<-X[,2]
  delta1<-delta[1]
  delta2<-delta[2]
  if(missing(eval.points)|any(is.na(eval.points))) {
    ax<-log(xlim+1/ngrid)
    ay<-log(ylim+1/ngrid)
    eval1<-exp(seq(ax[1],ax[2],length=ngrid))-1/ngrid
    eval2<-exp(seq(ax[1],ay[2],length=ngrid))-1/ngrid
   }
  else
   { eval1<-eval.points[,1];  eval2<-eval.points[,2]}
  pdf<- sm.density.eval.2d(log(x1+delta1),log(x2+delta2),
        h=h,hmult=1,xnew=log(eval1+delta1),ynew=log(eval2+delta2),
        eval.type=eval.type)
  est[is.na(est)]<-0
  if(eval.type=="points")
     est<-pdf$estimate/((eval1+delta1)*(eval2+delta2))
  else
     est<-pdf$estimate/outer(eval1+delta1,eval2+delta2)
  invisible(list(x1=eval1, x2=eval2, estimate=est,h=h))
}


sm.density.positive.grid<- function(X,
      h=c(hnorm(log(X[,1]+delta[1])),hnorm(log(X[,2]+delta[2]))),
      delta=c(min(X[,1]),min(X[,2])), ngrid=50, eval.points=NA,
      xlim=range(X[,1]),ylim=range(X[,2]))
{# converts outcome of sm.density.positive.2d to regular grid
   f <- sm.density.positive.2d(X,h,delta,ngrid=ngrid,
                   xlim=xlim,ylim=ylim,eval.type="grid")
   xx<-rep(f$x1,length(f$x2))
   yy<-rep(f$x2,rep(length(f$x2),length(f$x1)))
   zz<-as.vector(f$est,byrow=T)
   f.int<-interp(xx,yy,zz)
   invisible(list(eval.points = cbind(f.int$x, f.int$y),
   	estimate=f.int$z, h=h))
}


sm.density.3d <-
function(data, h=hnorm(data), hmult=1, contour = 75, ngrid = 20,
         shadow = T, xlab = NA, ylab = NA, zlab = NA, theta = pi/4, phi = pi/4,
         colour = T, title.colour = "magenta", label.colour = "magenta",
         axes.colour = "red",
         plot.colour = "blue", shadow.colour = "cyan", plt = T, cex = NA,
         maxpoly = 10000
	)
{
#
# ====================================================
# Function to plot a contour of a 3-d density estimate
# ====================================================
#
#	Written by Stuart Young.
#	Amended by Adrian Bowman.
#
# Input:
#   data	  = Data set. The function will also extract variable
#		    names, if the columns of the matrix are named.
#   contour	  = Contour level required.
#   gridsize	  = The dimension of the grid used in calculating the 3-d
#		    density estimate.  The smaller the gridsize, the
#		    "rougher" the final contour.
#   shadow	  = If True, a 'shadow' is drawn on the base of the
#		    perspective plot.
#   xlab	  = x-variable name (similarly for ylab and zlab)
#   theta	  = Horizontal angle of rotation for perspective plot.
#   phi	 	  = Vertical angle of rotation for perspective plot.
#   colour	  = If True, the plot is done in the various colours
#		    specified (see below).
#		    If False, only the figure itself is coloured (this over-
#		    rides any colours given in the parameters below).
#   title.colour  = Colour of the title (default is magenta).
#   label.colour  = Colour of the axis labels (default is magenta).
#   axes.colour   = Colour of the axes (default is red).
#   plot.colour	  = Colour of the contour (default is blue).
#   shadow.colour = Colour of the shadow (default is cyan).
#   plt		  = If True, then the plot is drawn.  If False, then the
#		    function simply returns the coordinate data for the
#		    3-d surface.  This over-rides other all other plot
#		    commands.
#   cex		  = Character expansion for the axis labels on the
#		    perspective plot.
#   maxpoly	  = This is required by the Fortran program for determining
#		    array sizes.  It sets the maximum number of polygons
#		    that can be stored.  The default value ought to be
#		    sufficient.  If the function crashes, increase the
#		    value of maxpoly.
#
#
# This function uses the Fortran subroutine 'npcontour'.  Before calling
# this function, the command dyn.load("fortran/npcontour.o") must be issued.
#
#
# Extract the three variables
#
	x <- data[, 1]
	y <- data[, 2]
	z <- data[, 3]	#
#
# Extract the variable names, as required
#
	if(is.na(xlab))
		if(!is.null(attributes(data)$dimnames)) {
			xlab <- attributes(data)$dimnames[[2]][1]
			ylab <- attributes(data)$dimnames[[2]][2]
			zlab <- attributes(data)$dimnames[[2]][3]
		}
		else {
			xlab <- ""
			ylab <- ""
			zlab <- ""
		}
#
# Calculate various summary statistics
#
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
	if(is.na(ngrid)) ngrid <- 12
	gridsize <- ngrid
#
# Calculate Normal optimal smoothing parameter
#
	hx <- h[1] * hmult
	hy <- h[2] * hmult
	hz <- h[3] * hmult
#
# Calculate density estimate at each data point
#
	fhat <- vector()
	for(i in (1:lng)) {
		fhat[i] <- sum(exp(-0.5 * (((x[i] - x)/hx)^2 + ((y[i] - y)/hy)^
			2 + ((z[i] - z)/hz)^2)))
	}
	fhat <- fhat/(lng * hx * hy * hz * (2 * pi)^(3/2))
	fhat <- sort(fhat)
	fmax <- fhat[lng]	#
#
# Find contour level (suitably scaled)
#
	p <- (lng * (100 - contour))/100
	p.trunc <- trunc(p)
	height <- fhat[p.trunc] + (p - p.trunc) * (fhat[p.trunc + 1] - fhat[
		p.trunc])
	height <- (height * 90)/fmax	#
#
# Calculate grid of x,y,z values
#
	xx <- seq(xmin, xmax, length = gridsize)
	yy <- seq(ymin, ymax, length = gridsize)
	zz <- seq(zmin, zmax, length = gridsize)	#
#
# Call Fortran subroutine to calculate coordinates
#
#	cat("\nCalculating surface coordinates\n")
	tmp <- as.double(rep(0, maxpoly * 3))
	result <- .Fortran("npcont",
		as.double(x),
		as.double(y),
		as.double(z),
		as.double(xx),
		as.double(yy),
		as.double(zz),
		as.double(hx),
		as.double(hy),
		as.double(hz),
		as.integer(lng),
		as.integer(gridsize),
		as.double(fmax),
		as.double(height),
		as.integer(maxpoly),
		as.double(rep(0, gridsize^3)),
		as.integer(0),
		as.double(rep(0, maxpoly * 3 * 3)),
		tmp,
		tmp,
		tmp)
	lng.coord <- 3 * result[[16]]
	xcoord <- result[[18]][1:lng.coord]
	ycoord <- result[[19]][1:lng.coord]
	zcoord <- result[[20]][1:lng.coord]
	coord <- cbind(xcoord, ycoord, zcoord)
	attributes(coord)$dimnames <- list(character(), c(xlab, ylab, zlab))	#
#
# Plot the contour, if requested
#
	if(plt) {
#		cat("\nPlotting will shortly commence\n\n")
		np.contour.plot.3d.(coord, data,shadow,
			gridsize, 3, xmin, xmax, ymin, ymax, zmin, zmax, xlab,
			ylab, zlab, theta, phi, colour, title.colour,
			label.colour, axes.colour, plot.colour, shadow.colour,
			cex)
	}
	else coord
}
"np.contour.plot.3d."<-
function(coord, data = matrix(), shadow = T,
	gridsize = 20, numpts = 3, xmin = NA, xmax = NA, ymin = NA, ymax = NA,
	zmin = NA, zmax = NA, xlab = NA, ylab = NA, zlab = NA, theta = pi/4,
	phi = pi/4, colour = F, title.colour = 3, label.colour = 3, axes.colour
	 = 6, plot.colour = "blue", shadow.colour = "cyan", cex = NA)
{
#
# ====================================================
# Function to plot a contour of a 3-d density estimate
# ====================================================
#
# Note: This function plots the 3-d coordinate data, which must be
#	supplied.
#
# Input:
#   coord	  = 3-column matrix of (x,y,z) coordinates of the polygons
#   data	  = Original data set (optional) - if the minimum and maximum
#		    values are not supplied, the function will extract them from
#		    the original data. The function will also extract variable
#		    names, if the columns of the matrix are named.
#   shadow	  = If True, a 'shadow' is drawn on the base of the perspective
#		    plot.
#   gridsize	  = The gridsize used to generate the polygons.
#   numpts	  = The number of vertices in each polygon.
#   xmin	  = Minimum observed value (from the data) for the x-variable.
#   xmax	  = Maximum observed value (similarly for ymin, ymax, zmin
#		    and zmax).
#   xlab	  = x-variable name (similarly for ylab and zlab).
#   theta	  = Horizontal angle of rotation for perspective plot.
#   phi		  = Vertical angle of rotation for perspective plot.
#   colour	  = If True, the plot is done in the various colours
#		    specified (see below).
#		    If False, only the figure itself is coloured (this over-
#		    rides any colours given in the parameters below).
#   title.colour  = Colour of the title (default is magenta).
#   label.colour  = Colour of the axis labels (default is magenta).
#   axes.colour   = Colour of the axes (default is red).
#   plot.colour	  = Colour of the contour (default is blue).
#   shadow.colour = Colour of the shadow (default is cyan).
#   cex		  = Character expansion for the axis labels on the
#		    perspective plot.
#
#
# Colours
#
#	1 = Yellow	2 = Cyan	3 = Magenta
#	4 = Green	5 = Blue	6 = Red
#
	col <- par("col")
	on.exit(par(col = col))
	if(!colour) {
		title.colour <- col
		label.colour <- col
		axes.colour <- col
	}
#
# Specify square plotting region
#
	opar <- par(pty = "s")
	on.exit(par(opar))	#
#
# Extract the variable names, as required
#
	if(is.na(xlab))
		if(!is.null(attributes(data)$dimnames)) {
			xlab <- attributes(data)$dimnames[[2]][1]
			ylab <- attributes(data)$dimnames[[2]][2]
			zlab <- attributes(data)$dimnames[[2]][3]
		}
		else {
			xlab <- ""
			ylab <- ""
			zlab <- ""
		}
#
# Extract min/max values and variable names, if not given
#
	if(is.na(xmin)) {
		xmin <- min(data[, 1])
		xmax <- max(data[, 1])
		ymin <- min(data[, 2])
		ymax <- max(data[, 2])
		zmin <- min(data[, 3])
		zmax <- max(data[, 3])
	}
#
# Store the scale matrix
#
	sc <- matrix(c(xmin, xmax, ymin, ymax, zmin, zmax), ncol = 2, byrow = T
		)	#
#
# Transform the coordinate data from the (1,gridsize) scale back to the
# original scale.
# Before scaling, the vector dst is calculated.  This records the distance
# of each point from the back of the axes system, and is used in the
# perspective plot below.
#
	dst <- sqrt(coord[, 1]^2 + coord[, 2]^2 + (zmax - coord[, 3])^2)
	coord[, 1] <- xmin + (coord[, 1] - 1)/(gridsize - 1) * (xmax - xmin)
	coord[, 2] <- ymin + (coord[, 2] - 1)/(gridsize - 1) * (ymax - ymin)
	coord[, 3] <- zmin + (coord[, 3] - 1)/(gridsize - 1) * (zmax - zmin)	#
#
# Create a matrix where each row represents one polygon, in the order
# (x1,y1,z1), (x2,y2,z2), ...
#
	poly <- matrix(t(coord), ncol = 3 * numpts, byrow = T)
	numpoly <- nrow(poly)	#
#
#
# ------------------------
# Plot 1: perspective plot
# ------------------------
#
# First, the polygons are ordered so that ones nearest the back of the
# axes system are plotted first.  The distance from each point to the back
# corner (1,1,gridsize) was calculated above and stored as the vector dst.
# Now, for each polygon, the minimum of these distances for that polygon
# is recorded.  The polygons are then ordered according to these values.
#
	dst <- matrix(dst, ncol = numpts, byrow = T)
	# Each row represents a polygon
	min.dst <- apply(dst, 1, min)	# calculate minimum for each polygon
	xyz <- poly[order(min.dst),  ]	# order by these minimum values
	par(col = axes.colour)
	box1.(theta, phi, sc, title.colour)
	if(shadow) {
		for(i in (1:numpoly)) {
			polyi <- matrix(xyz[i,  ], ncol = 3, byrow = T)
			polygon(coord2.(polyi[, 1], ymin, polyi[, 3], theta,
                                        phi, sc),
                                col = shadow.colour, border = shadow.colour)
	}
	}
#
# Draw the wire-frame figure.  Each polygon is drawn twice: the first time in
# background colour, to overwrite anything underneath (i.e. "behind" in a 3-d
# view).  The next is done with no shading, to create a wire-frame polygon.
#
	par(col = plot.colour)
	for(i in (1:numpoly)) {
		polyi <- matrix(xyz[i,  ], ncol = 3, byrow = T)
		polygon(coord2.(polyi[, 1], polyi[, 2], polyi[, 3], theta, phi,
			sc), col = 0)
		polygon(coord2.(polyi[, 1], polyi[, 2], polyi[, 3], theta, phi,
			sc))
	}
	par(col = axes.colour)
	box2.(theta, phi, sc, c(xlab, ylab, zlab), label.colour, cex)	#
#
#
}


"coord2."<-
function(x, y, z, theta, phi, sc)
{
#
# This function takes a set of 3-dimensional coordinates (x,y,z)
# {where (x,z)-axes are the base, and y is the upright axis}
# and returns a 2-dimensional matrix of corresponding (x',y')
# coordinates for plotting on a 2-d plane.
#
# Input : 	x,y,z 		- vectors of coordinates
#		theta, phi 	- angles of rotation
#		sc		- 2-column scale matrix, holding the values
#				  min(x), max(x)
#				  min(y), max(y)
#				  min(z), max(z)
#
	x <- -1 + ((x - sc[1, 1]) * 2)/(sc[1, 2] - sc[1, 1])
	y <- -1 + ((y - sc[2, 1]) * 2)/(sc[2, 2] - sc[2, 1])
	z <- -((-1 + ((z - sc[3, 1]) * 2)/(sc[3, 2] - sc[3, 1])))
	co <- cbind(x * cos(theta) - z * sin(theta), y * cos(phi) - (x * sin(
		theta) + z * cos(theta)) * sin(phi))
	co
}
"box1."<-
function(theta = pi/6, phi = pi/6, sc, col = par("col"),
	axes.lim = sc)
{
#
# Draw the back of the axes system
#
	plot(c( - sqrt(3), sqrt(3)), c( - sqrt(3), sqrt(3)), type = "n", xaxt
		 = "n", yaxt = "n", bty = "n", xlab = "", ylab = "")
	lines(rbind(coord2.(axes.lim[1, 1], axes.lim[2, 1], axes.lim[3, 1],
		theta, phi, sc), coord2.(axes.lim[1, 1], axes.lim[2, 1],
		axes.lim[3, 2], theta, phi, sc), coord2.(axes.lim[1, 2],
		axes.lim[2, 1], axes.lim[3, 2], theta, phi, sc)), lty = 3)
	lines(rbind(coord2.(axes.lim[1, 1], axes.lim[2, 1], axes.lim[3, 2],
		theta, phi, sc), coord2.(axes.lim[1, 1], axes.lim[2, 2],
		axes.lim[3, 2], theta, phi, sc)), lty = 3)
}


box2. <- function(theta = pi/6, phi = pi/6, sc, labels = c("", "", ""),
		col = par("col"), cex = 9/10, axes.lim = sc)
{
#
# Draw the front of the axes system
#
	lines(rbind(coord2.(axes.lim[1, 2], axes.lim[2, 1], axes.lim[3, 1],
		theta, phi, sc), coord2.(axes.lim[1, 2], axes.lim[2, 2],
		axes.lim[3, 1], theta, phi, sc), coord2.(axes.lim[1, 1],
		axes.lim[2, 2], axes.lim[3, 1], theta, phi, sc), coord2.(
		axes.lim[1, 1], axes.lim[2, 1], axes.lim[3, 1], theta, phi, sc),
		coord2.(axes.lim[1, 2], axes.lim[2, 1], axes.lim[3, 1], theta,
		phi, sc), coord2.(axes.lim[1, 2], axes.lim[2, 1], axes.lim[3, 2
		], theta, phi, sc), coord2.(axes.lim[1, 2], axes.lim[2, 2],
		axes.lim[3, 2], theta, phi, sc), coord2.(axes.lim[1, 1],
		axes.lim[2, 2], axes.lim[3, 2], theta, phi, sc), coord2.(
		axes.lim[1, 1], axes.lim[2, 2], axes.lim[3, 1], theta, phi, sc)
		), lty = 3)
	lines(rbind(coord2.(axes.lim[1, 2], axes.lim[2, 2], axes.lim[3, 1],
		theta, phi, sc), coord2.(axes.lim[1, 2], axes.lim[2, 2],
		axes.lim[3, 2], theta, phi, sc)), lty = 3)	#
#
# Label the axes
#
	cex <- 0.9
	text(rbind(coord2.((axes.lim[1, 1] + axes.lim[1, 2])/2, axes.lim[2, 1],
		axes.lim[3, 1], theta, phi, sc)), labels[1], adj = 1, cex = cex,
		col = col)
	text(rbind(coord2.(axes.lim[1, 1], (axes.lim[2, 1] + axes.lim[2, 2])/2,
		axes.lim[3, 1], theta, phi, sc)), labels[2], adj = 1, cex = cex,
		col = col)
	text(rbind(coord2.(axes.lim[1, 2], axes.lim[2, 1], (axes.lim[3, 1] +
		axes.lim[3, 2])/2, theta, phi, sc)), labels[3], adj = 0, cex =
		cex, col = col)
}
sm.density.compare <- function(x, group, h = NA, display = "lines",
			model = "none", test = T, band = T, ngrid = 50,
			nboot = 100, monitor = T,
			xlab, ylab, xlim, lty, ...)

{

y <- x

if (model=="none") {
	band <- F
	test <- F
	}
if (display == "none") band <- F

fact        <- factor(group)
fact.levels <- levels(fact)
nlevels     <- length(fact.levels)
ni <- table(fact)

if (band & (nlevels > 2)) {
  cat("Reference band available to compare two groups only.","\n")
  band <- F
  }

if (missing(xlab)) xlab <- deparse(substitute(x))
if (missing(ylab)) ylab <- "Density"
if (missing(xlim)) xlim <- c(min(y) - diff(range(y))/4,
			max(y) + diff(range(y))/4)
if (is.na(h)) {
    hvec <- tapply(y, fact, FUN = "hnorm")
    h <- exp(mean(log(hvec)))
    }

if (missing(lty)) lty <- 1:nlevels

estimate <- matrix(0, ncol = ngrid, nrow = nlevels)
se       <- matrix(0, ncol = ngrid, nrow = nlevels)
for (i in 1:nlevels) {
  sm <- sm.density(y[fact == fact.levels[i]], h = h,
  			xlim = xlim, ngrid = ngrid, display = "none")
  estimate[i,] <- sm$estimate
  se[i,]       <- sm$se
  }
eval.points <- sm$eval.points

if (!(display == "none" | band)) {
  plot(xlim, c(0, 1.1 * max(as.vector(estimate))),
		xlab = xlab, ylab = ylab, type = "n", ...)
  for (i in 1:nlevels) lines(eval.points, estimate[i,], lty = lty[i])
  }

est <- NULL
p   <- NULL
if (model == "equal" & test) {

if (nlevels == 2) {ts <- sum((estimate[1,] - estimate[2,])^2)}
else {
  sm.mean <- sm.density(y, h = h, xlim = xlim, ngrid = ngrid,
			display = "none")$estimate
  ts <- 0
  for (i in 1:nlevels) ts <- ts + ni[i]*sum((estimate[i,] - sm.mean)^2)
  }

p <- 0
est.star <- matrix(0, ncol = ngrid, nrow = nlevels)
for (i in 1:nboot) {

  ind <- (1:length(y))
  for (i in 1:nlevels) {
    indi <- sample((1:length(ind)), ni[i])
    est.star[i,] <- sm.density(y[ind[indi]], h = h,
		ngrid = ngrid, xlim = xlim, display = "none")$estimate
    ind <- ind[-indi]
    }

  if (nlevels == 2) {ts.star <- sum((est.star[1,] - est.star[2,])^2)}
  else {
    sm.mean <- sm.density(y, h = h, xlim = xlim, ngrid = ngrid,
				display = "none")$estimate
    ts.star <- 0
    for (i in 1:nlevels) {ts.star <- ts.star + ni[i]*sum((est.star[i,] - sm.mean)^2)}
    }

  if (ts.star > ts) p <- p + 1
  if (monitor) {
    cat(i)
    cat(" ")
    }
  }
p <- p/nboot
cat("\n", "Test of equal densities:  p-value = ", p, "\n")
est <- list(p = p, h = h)
}

if (model == "equal" & band ) {
  av <- (sqrt(estimate[1,]) + sqrt(estimate[2,]))/2
  se <- sqrt(se[1,]^2 + se[2,]^2)
  upper <- (av + se)^2
  lower <- pmax(av - se, 0)^2
  plot(xlim, c(0, 1.1 * max(as.vector(estimate), upper)),
		xlab = xlab, ylab = ylab, type = "n", ...)
  polygon(c(eval.points, rev(eval.points)), c(upper, rev(lower)), col = "cyan",
		border = F)
  lines(eval.points, estimate[1,], lty = lty[1])
  lines(eval.points, estimate[2,], lty = lty[2])
  est <- list(p = p, upper = upper, lower = lower, h = h)
  }

invisible(est)

}
# version 17-12-1996

# parameter: display=none|estimate|se, add=T|F

sm.logit <- function(x, y, N=rep(1,length(y)), h, ngrid=25,
        eval.points, add=F, display="estimate",
        xlab, ylab, pch=1, col=2, ...)
{
  xlab0 <- deparse(substitute(x)); xlab0
  ylab0 <- deparse(substitute(y)); ylab0
  logistic<-function(x) 1/(1+exp(-x))
  logit <-function(p) log(p/(1-p))
  y <- as.integer(y)
  if(min(diff(x))<0) {y <- y[order(x)]; x <- sort(x)}
  n <- length(y)
  if(missing(eval.points)) eval.points<-seq(min(x),max(x),length=ngrid)
  if(all(N==1)) yplot<-jitter(y, amount=0) else yplot<-y/N
  if(display!="none" & add==F) {
    if(missing(xlab)) xlab<- xlab0
    if(missing(ylab)) ylab<- paste("Pr{",ylab0,"}",sep="")
    plot(x,yplot,ylim=c(0,1), xlab=xlab, ylab=ylab, ..., col=1, type="n")
    abline(0,0, col=1, lty=3)
    abline(1,0, col=1, lty=3)
    }
  if(display!="none")  points(x,yplot,pch=pch,col=col)
  neval <- length(eval.points)
  # W <- sm.weight(x,eval.points,h,poly.index=0)
  W <- matrix(rep(eval.points, rep(n, neval)), ncol = n, byrow = T)
  W <- W - matrix(rep(x, neval), ncol = n, byrow = T)
  W <- exp(-0.5 * (W/h)^2)
  cat("Cycles per point: ")
  est <- LP <- st.err <- rep(NA,neval)
  Y <- cbind(c(y,0),c(N,1)-c(y,0))
  X <- cbind(rep(1,n+1),c(x,0))
  start<-c(log((y+0.5)/(N-y+1)),0)
  for(k in 1:neval){
    X[n+1,2]<-eval.points[k]
    dimnames(X) <- list(NULL, c("x", "y"))
    fit<-glm.fit(X, Y, c(W[k,],0), family=binomial(), etastart=start)
    start<-fit$linear.predictors
    LP[k]<-start[n+1]
    st.err[k] <-
       sqrt(X[n+1,] %*% solve(crossprod(fit$R)) %*% as.matrix(X[n+1,]))
    cat(fit$iter);cat(" ")
  } #end for
  cat("\n")
  est<-logistic(LP)
  result<-list(eval.points=eval.points,estimate=est,linear.predictors=LP)
  result$upper<-logistic(LP+2*st.err)
  result$lower<-logistic(LP-2*st.err)
  result$se<-st.err
  if(display!="none"){
    lines(result$eval.points,result$estimate,col=col)
    if(display=="se") {
      lines(result$eval.points,result$lower,lty=3,col=col)
      lines(result$eval.points,result$upper,lty=3,col=col)
      }
    }
  invisible(result)
}# end def fn



sm.logit.bootstrap <- function(x, y, N=rep(1,length(x)), h, nboot=100,
          degree=1, ...)
{
  logistic<-function(x) 1/(1+exp(-x))
  D <- function(y,N,mu) {
          p<-y/N
          d <- y*log(p*(1-mu)/((1-p)*mu))+N*log((1-p)/(1-mu))
          d[y==0] <- (-N*log(1-mu[y==0]))
          d[y==N] <- (-N*log(mu[y==N]))
          dev<-2*sum(d)
          if(is.na(dev)) dev<-Inf
          dev
          }
  n <- length(x)
  sm <- sm.logit(x,y,N,h,xlab=deparse(substitute(x)),
           ylab= paste("Pr{",deparse(substitute(y)),"}",sep=""), ...)
  X <- as.matrix(cbind(1,poly(x,degree)))
  dimnames(X) <- list(NULL, paste("V",1:ncol(X), sep=""))
  glm.model <- glm.fit(X, cbind(y,N-y), family=binomial())
  glm.fitted<- fitted(glm.model)
  lines(x,glm.fitted,lty=2,col=2)
  p.boot <- 0
  sm.orig <- sm.logit(x, y, N, h, eval.points=x, display="none")
  sm.fitted <- sm.orig$estimate
  #disp <- D(y,N,sm.fitted)/(n-degree-1)
  disp <- 1
  ts.orig <- (D(y,N,glm.fitted)-D(y,N,sm.fitted))/disp
  type("Dispersion parameter = ",disp)
  type("Test statistic = ", ts.orig)
  disp.par<-rep(NA,nboot)
  for (i in 1:nboot) {
    yboot<-rbinom(length(glm.fitted),N,glm.fitted)
    sm.fitted <- sm.logit(x,yboot,N,h,eval.points=x,display="none")$estimate
    # disp <- D(yboot,N,sm.fitted)/(n-degree-1)
    disp<-1
    if(disp<Inf) ts.boot <- (D(yboot,N,glm.fitted)-D(yboot,N,sm.fitted))/disp
            else {ts.boot<- -Inf ; cat("Inf deviance\n")}
    if (ts.boot > ts.orig) p.boot <- p.boot + 1
    lines(x, sm.fitted, lty=2,col=6)
    print(paste(i,p.boot,ts.boot,ts.orig))
    disp.par[i]<-disp
    }
  lines(sm$eval.points,sm$estimate)
  p.boot <- p.boot/nboot
  type("Observed significance = ", p.boot)
  invisible(list(test.statistic=ts.orig, significance=p.boot, disp=disp.par))
}


nnbr <- function(x, k) {

	 if (isMatrix(x)) {
	 	ndim <- 2
	 	n <- nrow(x)
	 	}
	 else {
	 	ndim <- 1
	 	n <- length(x)
	 	}

         knn <- vector("numeric", n)

         if (ndim==1) {
	    for (i in 1:length(x)) knn[i] <- sort(abs(x-x[i]))[k+1]
	    }
         if (ndim==2) {
            for(i in 1:length(x[, 1]))
                knn[i] <- sort(sqrt(((x[, 1] - x[i, 1])^2)/var(x[, 1]) +
                	((x[,2] - x[i, 2])^2)/var(x[, 2])))[k + 1]
            }

         knn
	 }

# version 30-1-1997


sm.poisson <- function(x, y,  h, ngrid=25, eval.points, add=F,
              display="estimate", xlab, ylab, pch=1, col=2, ...)
{
  xlab0 <- deparse(substitute(x)); xlab0
  ylab0 <- deparse(substitute(y)); ylab0
  y <- as.integer(y)
  if(min(diff(x))<0) {y <- y[order(x)]; x <- sort(x)}
  n <- length(y)
  if(missing(eval.points)) eval.points<-seq(min(x),max(x),length=ngrid)
  if(display!="none" & add==F) {
    if(missing(xlab)) xlab<- ylab0
    if(missing(ylab)) ylab<- ylab0
    plot(x, y, xlab=xlab, ylab=ylab, ..., col=1, type="n")
    }
  if(display!="none")  points(x,y,pch=pch,col=col)
  neval <- length(eval.points)
  W <- matrix(rep(eval.points, rep(n, neval)), ncol = n, byrow = T)
  W <- W - matrix(rep(x, neval), ncol = n, byrow = T)
  W <- exp(-0.5 * (W/h)^2)
  est <- rep(NA,neval)
  cat("Cycles per point: ")
  LP<-st.err<-rep(NA,neval)
  X <- cbind(v1=rep(1,n+1),v2=c(x,0))
  Y <- c(y,0)
  start<-c(log(y+0.5),0)
  for(k in 1:neval){
    X[n+1,2]<-eval.points[k]
    fit<-glm.fit(X, Y, w=c(W[k,],0), family=poisson(), etastart=start, ...)
    start<-fit$linear.predictors
    LP[k]<-start[n+1]
    st.err[k] <- sqrt(X[n+1,] %*% solve(crossprod(fit$R)) %*% as.matrix(X[n+1,]))
    cat(fit$iter);cat(" ")
  } #end for
  cat("\n")
  est<-exp(LP)
  result<-list(eval.points=eval.points,estimate=est,linear.predictors=LP)
  result$upper<-exp(LP+2*st.err)
  result$lower<-exp(LP-2*st.err)
  result$st.err<-st.err
  if(display!="none"){
    lines(result$eval.points,result$estimate,col=col)
    if(display=="se") {
      lines(result$eval.points,result$lower,lty=3,col=col)
      lines(result$eval.points,result$upper,lty=3,col=col)
      }
    }
  invisible(result)
}# end def fn



sm.poisson.bootstrap <- function(x, y, h, nboot=100, degree=1, ...)
{
  D<-function(y,mu){
       d<-(mu-y+y*log(y/mu))
       d[y==0]<-mu[y==0]
       return(2*sum(d))
       }
  xlab0 <- deparse(substitute(x)); xlab0
  ylab0 <- deparse(substitute(y)); ylab0
  y <- as.integer(y)
  y <- y[order(x)]
  x <- sort(x)
  X <- cbind(v1=1,v2=poly(x,degree))
  sm <- sm.poisson(x, y, h, xlab=xlab0, ylab=ylab0, ...)
  glm.model <- glm.fit(X, y, family=poisson())
  glm.fitted<- fitted(glm.model)
  lines(x,glm.fitted,lty=2,col=2)
  p.boot <- 0
  sm.orig <- sm.poisson(x, y, h, eval.points=x, display="none")
  sm.fitted <- sm.orig$estimate
  disp <- D(y,sm.fitted)/(length(y)-degree-1)
  ts.orig <- (D(y,glm.fitted)-D(y,sm.fitted))/1  # disp
  type("Test statistic = ", ts.orig)
  disp.par<-rep(NA,nboot)
  for (i in 1:nboot) {
    yboot<-rpois(length(glm.fitted),glm.fitted)
    sm.fitted <- sm.poisson(x, yboot, h, eval.points=x, display="none")$estimate
    disp <- D(yboot,sm.fitted)/(length(y)-degree-1)
    ts.boot <- (D(yboot,glm.fitted)-D(yboot,sm.fitted))/1   #  disp
    if (ts.boot > ts.orig) p.boot <- p.boot + 1
    lines(x, sm.fitted, lty=2, col=6)
    print(paste(i,p.boot,ts.boot,ts.orig))
    disp.par[i]<-disp
    }
  lines(sm$eval.points,sm$estimate)
  p.boot <- p.boot/nboot
  type("Observed significance = ", p.boot)
  invisible(list(test.statistic=ts.orig, significance=p.boot, disp=disp.par))
}


#		Calculation of the probability that a quadratic
#		form y^TAY is greater than zero, where y~N(0,Sigma)

p.quad.moment <- function(A, Sigma) {
   B  <- A %*% Sigma
   k1 <- sum(diag(B))
   C <- B %*% B
   k2 <- 2 * sum(diag(C))
   k3 <- 8 * sum(diag(C %*% B))
   aa <- abs(k3/(4 * k2))
   bb <- (8 * k2^3)/k3^2
   cc <- k1 - aa * bb
   1 - pchisq( - cc/aa, bb)
   }
#     Nonparametric regression

sm.regression <- function(x, y, h, design.mat = NA, hmult = 1,
	h.weights = NA, poly.index = 1,
	model = "none", band = T, test = T, display = "lines", add = F,
	ngrid = NA, eval.points = NA, weights=rep(1,n),
	xlab = NA, ylab = NA, zlab = NA, hull = T, panel=F, lty=1, col=1,
	eye.mult = c(-6, -8, 5), ... ) {

if (length(dim(x))>0)
     { ndim <- 2 ; n<-dim(x)[1]}
  else
     { ndim <- 1; n<-length(x)}

if (ndim==1) {
	if (is.na(xlab))  xlab <- deparse(substitute(x))
	if (is.na(ylab))  ylab <- deparse(substitute(y))
	if (is.na(ngrid)) ngrid <- 50
	if (any(is.na(h.weights))) h.weights <- rep(1, length(x))
  	est <- sm.regression.1d(x, y, h, design.mat, hmult, h.weights,
  		  poly.index,
	          model, band, test, display, add, ngrid, eval.points,
	          weights, xlab, ylab, panel, lty, col, ... )
		}
else {
	if (is.na(ngrid)) ngrid <- 20
	if (is.na(xlab)) {
		if(!is.null(attributes(x)$dimnames))
		   xlab <- attributes(x)$dimnames[[2]][1]
		else xlab <- paste(deparse(substitute(x)),"[1]")
		}
	if (is.na(ylab)) {
		if(!is.null(attributes(x)$dimnames))
		   ylab <- attributes(x)$dimnames[[2]][2]
		else ylab <- paste(deparse(substitute(x)),"[2]")
		}
	if (is.na(zlab)) zlab <- deparse(substitute(y))
	if (any(is.na(h.weights))) h.weights <- rep(1, nrow(x))
	est <- sm.regression.2d(x, y, h, hmult, h.weights, poly.index,
	          model, test, display, ngrid, eval.points, weights,
	          xlab, ylab, zlab, hull, eye.mult, ... )
	}

invisible(est)

}

sm.regression.1d <- function(x, y, h, design.mat = NA,
		hmult = 1, h.weights = rep(1,length(x)), poly.index = 1,
	        model = "none", band = T, test = T, display = "lines", add = F,
	        ngrid = 50, eval.points = NA, weights = rep(1,length(x)),
	        xlab, ylab, panel, lty=1, col=1, ... )  {

if (model=="none") {
	band <- F
	test <- F
	}
if (add | display=="none") panel <- F
if (!(model=="none") & panel == F) test <- T

r <- list(x=NA, y=NA, model.y=NA, se=NA, sigma=NA,
		h=h*hmult, hweights=h.weights)
if (!add & !(display=="none")) plot(x, y, xlab=xlab, ylab=ylab, ...,
		lty=1, col=1)
if (!(display=="none")) r <- plot.regression(x,y,design.mat,h,hmult,h.weights,r,
				model,band,test=F,display,poly.index,
				ngrid, add = add, weights, lty, col)
if (test) rtest <- sm.regression.test(x, y, design.mat, h, hmult, h.weights,
			poly.index, model)


if (panel) {
	items <-      c("Bandwidth:",
#			"  - Cross-validation",
#			"  - Plug-in",
			"  - increase",
			"  - decrease",
			"  - movie up",
			"  - movie down",
#                        "Add linear band",
			"Exit")
#	if (band) items[8] <- "Remove linear band"

	ind <- menu(items, graphics=T, title="Nonparametric regression")

	while (items[ind]!="Exit") {
		if (items[ind]=="  - increase") {
			hmult <- hmult * 1.1
			}
		else if (items[ind]=="  - decrease") {
			hmult <- hmult / 1.1
			}
		else if (items[ind]=="  - movie up") {
			for (i in 1:6) {
				hmult <- hmult * 1.1
				r <- plot.regression(x,y,design.mat,
						h,hmult,h.weights,r,
						model,band, test=F, display,
						poly.index, ngrid,
					add = add, weights = weights, lty, col)
				}
			hmult <- hmult * 1.1
			}
		else if (items[ind]=="  - movie down") {
			for (i in 1:6) {
				hmult <- hmult / 1.1
				r <- plot.regression(x,y,design.mat,
						h,hmult,h.weights,r,
						model,band, test=F, display,
						poly.index,ngrid,
						add = add,
						weights = weights)
				}
			hmult <- hmult / 1.1
			}
		else if (items[ind]=="Add linear band" |
				items[ind]=="Remove linear band") {
			bandflag <- !bandflag
			if (!bandflag) polygon(c(r$x,rev(r$x)),
                   		c(mean(y)-2*r$se,mean(y)+2*rev(r$se)),col=0)
			if (items[ind]=="Add linear band") {
				items[ind] <- "Remove linear band"  }
			else (items[ind] <- "Add linear band")
			}

		r <- plot.regression(x,y,design.mat,h,hmult,h.weights,r,
				model,band, test=F, display, poly.index, ngrid,
				add = add, weights = weights, lty, col)
		cat("h = ", signif(h*hmult,7), "\n")

		ind <- menu(items, graphics=T, title="Nonparametric regression")
		}
	}

if (!(any(is.na(eval.points)))) r <- sm.regression.eval.1d(x,y,design.mat,
			h,hmult,h.weights,model,band,test,poly.index,
			ngrid,eval.points, weights)
   else if ((display=="none") & (model=="none")) {
		r <- sm.regression.eval.1d(x,y,design.mat,
			h,hmult,h.weights,model,band,test,poly.index,
			ngrid, eval.points=seq(min(x),max(x),length=ngrid),
			weights)
			}

if (test) r <- list(eval.points=r$eval.points, estimate=r$estimate,
   		model.y=r$model.y, se=r$se, sigma=r$sigma, h=r$h,
   		hweights=r$hweights, model=rtest$model, p=rtest$p)

r

}

#----------------------------Display graphics----------------------------

plot.regression <- function(x,y,design.mat,h,hmult,h.weights,r,
				model,band,test,display,poly.index=1,
				ngrid, add, weights, lty, col) {

  rnew <- sm.regression.eval.1d(x,y,design.mat,h,hmult,h.weights,
  			model,band,test,poly.index, ngrid,
  			weights = weights)
  if (length(r$eval.points) > 0 && !any(is.na(r$eval.points))) {
  	if (band) {
  		upper <- r$model.y + 2*r$se
  		upper <- pmin(pmax(upper,par()$usr[3]),par()$usr[4])
  		lower <- r$model.y - 2*r$se
  		lower <- pmin(pmax(lower,par()$usr[3]),par()$usr[4])
  		polygon(c(r$eval.points,rev(r$eval.points)),
  			c(lower,rev(upper)), col=0, border=F)
  		}
  	if (display == "se") {
  		upper <- r$estimate + 2*r$se
  		upper <- pmin(pmax(upper,par()$usr[3]),par()$usr[4])
  		lower <- r$estimate - 2*r$se
  		lower <- pmin(pmax(lower,par()$usr[3]),par()$usr[4])
  		lines(r$eval.points, upper, lty = 3, col = 0)
  		lines(r$eval.points, lower, lty = 3, col = 0)
  		}
  	lines(r$eval.points, r$estimate, col=0)
  	}
  if (band) {
  	upper <- rnew$model.y + 2*rnew$se
  	upper <- pmin(pmax(upper,par()$usr[3]),par()$usr[4])
  	lower <- rnew$model.y - 2*rnew$se
  	lower <- pmin(pmax(lower,par()$usr[3]),par()$usr[4])
  	polygon(c(rnew$eval.points,rev(rnew$eval.points)),
  		c(lower,rev(upper)), col="cyan", border=F)
  	}
  lines(rnew$eval.points, rnew$estimate, lty=lty, col=col)
  if ((model == "none") & (display == "se")) {
  	upper <- rnew$estimate + 2*rnew$se
  	upper <- pmin(pmax(upper,par()$usr[3]),par()$usr[4])
  	lower <- rnew$estimate - 2*rnew$se
  	lower <- pmin(pmax(lower,par()$usr[3]),par()$usr[4])
  	lines(rnew$eval.points, upper, lty = 3, col=col)
  	lines(rnew$eval.points, lower, lty = 3, col=col)
	}
  if (!add) points(x, y,col=1)
  box(col=1, lty=1)
  rnew
  }


#------------------regression estimate---------------------------------

sm.regression.eval.1d <- function(x, y, design.mat, h,
		 hmult=1, h.weights = rep(1,length(x)),
		 model="none", band=F, test=F,
                 poly.index = 1,
                 ngrid=50, eval.points=seq(min(x),max(x),length=ngrid),
                 weights = rep(1,length(x)))
  {
	w   <- sm.weight(x, eval.points, h, hmult, h.weights,
			poly.index, weights = weights)
        est <- as.vector(w %*% y)
        sig <- sm.sigma(x,y)

        n  <- length(x)
        ne <- length(eval.points)

        if (model=="none") {
        	model.y <- est
        	se  <- as.vector(sig * sqrt(((w^2) %*% rep(1,n))))
        	}
        else if ((model=="no.effect") | (model=="no effect")) {
        	if (is.na(as.vector(design.mat)[1])) {
        	     X <- matrix(rep(1,n),ncol=1)
        	     model.y <- rep(mean(y),ne)
        	     }
        	else {
        	     X <- design.mat
        	     model.y <- rep(0,ne)
        	     }
        	X <- diag(n) - X %*% solve(t(X) %*% X) %*% t(X)
        	se <- sig * sqrt(diag(w %*% X %*% t(w)))
        	}
        else if (model=="linear") {
        	e <- cbind(rep(1,ne), eval.points - mean(x))
        	l <- cbind(rep(1,n ), x - mean(x))
        	l <- e %*% solve(t(l) %*% l) %*% t(l)
        	model.y <- as.vector(l %*% y)
        	se <- as.vector(sig * sqrt(((w-l)^2) %*% rep(1,n)))
		}

        list(eval.points=eval.points, estimate=est, model.y=model.y, se=se, sigma=sig,
        	h=h*hmult, hweights=h.weights, weights=weights)

        }

#----------------------------Tests-----------------------------------

sm.regression.test <- function(x, y, design.mat = NA, h,
			hmult = 1, h.weights = rep(1,length(x)), poly.index=1,
			model = "no.effect")
  {
	if (length(dim(x))>0) {
     		ndim <- 2
     		n<-dim(x)[1]
		W <- sm.weight2(x, x, h, hmult, h.weights, poly.index)
		S <- cbind(rep(1,n), x[,1] - mean(x[,1]), x[,2] - mean(x[,2]))
     		}
  	else {
     		ndim <- 1
     		n<-length(x)
        	n <- length(x)
		W <- sm.weight(x, x, h, hmult, h.weights, poly.index)
		S <- cbind(rep(1,n), x - mean(x))
		}

        if ((model=="no.effect") | (model=="no effect")) {
        	if (is.na(as.vector(design.mat)[1]))
        	   S <- matrix(rep(1,n),ncol=1)
        	   else S <- design.mat
        	S  <- diag(n) - S %*% solve(t(S) %*% S) %*% t(S)
		W  <- diag(n) - W
		W  <- t(W) %*% W
        	e  <- as.vector(S %*% y)
        	r0 <- sum(e^2)
        	r1 <- as.numeric(t(e) %*% W %*% e)
        	ts <- (r0-r1)/r1
        	p  <- p.quad.moment(diag(n)-(1+ts)*W, S)
        	}

        else if (model=="linear") {
        	S  <- diag(n) - S %*% solve(t(S) %*% S) %*% t(S)
		W  <- diag(n) - W
		W  <- t(W) %*% W
        	e  <- as.vector(S %*% y)
        	r0 <- sum(e^2)
        	r1 <- as.numeric(t(e) %*% W %*% e)
        	ts <- (r0-r1)/r1
        	p  <- p.quad.moment(diag(n)-(1+ts)*W, S)
        	}

	print(paste("Test of", model,"model:  significance = ",round(p,3)))
        list(model = model, p = p, h = h*hmult, hweights = h.weights)

        }

#-------------------Estimation of error standard deviation---------------

sm.sigma <- function(x, y, diff.ord = 2) {

  if (diff.ord == 1) {
    yy <- y[order(x)]
    sig <- sqrt(sum((yy[2:length(yy)]-yy[1:(length(yy)-1)])^2)/(2*(length(yy)-1)))
    }

  else {
    n   <- length(x)
    yy  <- y[order(x)]
    xx  <- sort(x)
    xx1 <- diff(xx)
    xx2 <- diff(xx, lag = 2)

    a   <- xx1[-1]/xx2
    b   <- xx1[-(n-1)]/xx2
    a[xx2==0] <- 0.5
    b[xx2==0] <- 0.5
    cc  <- sqrt(a^2 + b^2 + 1)
    eps <- yy[1:(n-2)] * a/cc + yy[3:n] * b/cc - yy[2:(n-1)]/cc

    sig <- sqrt(sum(eps^2)/(n-2))
    }

  sig

  }


#-----------------Weight matrix for nonparametric regression---------------

sm.weight <- function(x, eval.points, h, hmult=1, h.weights = rep(1,length(x)),
		poly.index=1, cross = F, weights = rep(1,length(x))) {

n  <- length(x)
ne <- length(eval.points)

if (poly.index==0) {
        wd  <- matrix(rep(eval.points, rep(n, ne)), ncol = n, byrow = T)
        wd  <- wd - matrix(rep(x, ne), ncol = n, byrow = T)
        w   <- matrix(rep(h.weights,ne), ncol = n, byrow = T)
        w   <- exp(-.5 * (wd/(h*hmult*w))^2)
        w   <- w * matrix(rep(weights,ne), ncol = n, byrow = T)
        if (cross) diag(w) <- 0
        den <- w %*% rep(1, n)
        w   <- w / matrix(rep(den,n), ncol=n)
        }
else if (poly.index==1) {
	wd <- matrix(rep(eval.points, rep(n, ne)), ncol = n, byrow = T)
        wd <- wd - matrix(rep(x, ne), ncol = n, byrow = T)
        w  <- matrix(rep(h.weights,ne), ncol = n, byrow = T)
        w  <- exp(-.5 * (wd/(h*hmult*w))^2)
        w  <- w * matrix(rep(weights,ne), ncol = n, byrow = T)
        if (cross) diag(w) <- 0

	s0 <-  w         %*% rep(1,n)
	s1 <- (w * wd)   %*% rep(1,n)
	s2 <- (w * wd^2) %*% rep(1,n)

	w  <-  w * (matrix(rep(s2,n),ncol=n) - wd * matrix(rep(s1,n),ncol=n))
	w  <-  w / (matrix(rep(s2,n),ncol=n) * matrix(rep(s0,n),ncol=n)
			- matrix(rep(s1,n),ncol=n)^2)
	}

}

sm.regression.2d <- function(x, y, h, hmult, h.weights = NA,
		poly.index = 1, model = "none", test = T, display = "persp",
		ngrid = 20, eval.points = NA, weights = rep(1,n),
		xlab = NA, ylab = NA, zlab = NA, hull = T,
		eye.mult = c(-6, -8, 5), ...)

#	Evaluate a 2-d nonparametric regression estimate over a square grid

{

if (model=="none") test <- F

n <- length(y)
if (any(is.na(h.weights))) h.weights <- rep(1,n)

ev.points <- NA
est       <- NA

x1grid      <- seq(min(x[,1]), max(x[,1]), length=ngrid)
x2grid      <- seq(min(x[,2]), max(x[,2]), length=ngrid)
e.points <- cbind(x1grid, x2grid)

if (!(display == "none")) {
   est <- sm.regression.eval.2d(x, y, h, hmult, h.weights, model,
   		poly.index, ngrid, e.points, hull)
   ev.points <- e.points
   eye <- c(x1grid[1] + eye.mult[1] * diff(range(x1grid)),
       		x2grid[1] + eye.mult[2] * diff(range(x2grid)),
       		max(est[!is.na(est)])  +
       			eye.mult[3] * diff(range(est[!is.na(est)])))
#   persp(x1grid, x2grid, est, xlab = xlab, ylab = ylab, zlab = zlab,
#       		eye = eye, ...)
   persp(x1grid, x2grid, est, theta = -30, phi = 40, d = 4)
   }

if (!(any(is.na(as.vector(eval.points))))) {
   eval.type <- "points"
   ev.points <- eval.points
   w   <- sm.weight2(x, eval.points, h, hmult, h.weights, poly.index)
   est <- matrix(w %*% y, ncol = ngrid)
   }
   else if ((display=="none") & (model=="none")) {
   	   est <- sm.regression.eval.2d(x, y, h, hmult, h.weights, model,
   		poly.index, ngrid, e.points, hull)
   	   ev.points <- e.points
	   }


model.y <- NA
se      <- NA
sigma   <- NA

r <- list(eval.points = ev.points, estimate = est, model.y = model.y,
	  se = se, sigma = sigma, h = h * hmult, hweights = h.weights)

if (test) {
	rtest <- sm.regression.test(x, y, design.mat=NA, h, hmult,
			h.weights, poly.index, model)
	r <- list(eval.points = ev.points, estimate = est,
		model.y = model.y, se = se, sigma = sigma, h = h * hmult,
		hweights = h.weights, model=rtest$model, p=rtest$p)
	}

r

}


#-----------------Regression function evaluated over a grid----------------

sm.regression.eval.2d <- function(x, y, h, hmult, h.weights, model,
   		poly.index, ngrid, eval.points, hull = T)

{

n     <- nrow(x)
ngrid <- nrow(eval.points)

if (any(is.na(h.weights))) h.weights <- rep(1,n)

wd1   <- matrix(rep(eval.points[,1], n), ncol = n)
wd1   <- wd1 - matrix(rep(x[,1], ngrid), ncol = n, byrow = T)
wd2   <- matrix(rep(eval.points[,2], n), ncol = n)
wd2   <- wd2 - matrix(rep(x[,2], ngrid), ncol = n, byrow = T)
wy    <- matrix(rep(h.weights, ngrid), ncol = n, byrow = T)
w1    <- exp(-.5 * (wd1 / (h[1] * hmult * wy))^2)
w2    <- exp(-.5 * (wd2 / (h[2] * hmult * wy))^2)
wy    <- matrix(rep(y, ngrid), ncol = n, byrow=T)

if (poly.index == 0) {

   est   <- w1 %*% t(w2 * wy) / (w1 %*% t(w2))

   }

if (poly.index == 1) {

   a11   <-  w1          %*% t(w2)
   a12   <- (w1 * wd1)   %*% t(w2)
   a13   <-  w1          %*% t(w2 * wd2)
   a22   <- (w1 * wd1^2) %*% t(w2)
   a23   <- (w1 * wd1)   %*% t(w2 * wd2)
   a33   <-  w1          %*% t(w2 * wd2^2)

   d     <- a22 * a33 - a23^2

   b1    <- 1 / (a11 - ((a12*a33 - a13*a23)*a12 + (a13*a22 - a12*a23)*a13)/d)
   b2    <- (a13*a23 - a12*a33) * b1 / d
   b3    <- (a12*a23 - a13*a22) * b1 / d

   c1    <-  w1        %*% t(w2 * wy)
   c2    <- (w1 * wd1) %*% t(w2 * wy)
   c3    <-  w1        %*% t(w2 * wy * wd2)

   est   <- b1 * c1 + b2 * c2 + b3 * c3

   }

if (hull) {
  hull.points <- x[order(x[,1], x[,2]),]
  dh          <- diff(hull.points)
  hull.points <- hull.points[c(T, !((dh[,1]==0) & (dh[,2]==0))),]
  hull.points <- hull.points[chull(hull.points),]
  nh          <- nrow(hull.points)
  gstep       <- matrix(rep(eval.points[2,] - eval.points[1,], nh),
  			ncol = 2, byrow = T)
  hp.start    <- matrix(rep(eval.points[1,], nh), ncol = 2, byrow = T)
  hull.points <- hp.start +
  		gstep * round((hull.points - hp.start) / gstep)
  hull.points <- hull.points[chull(hull.points),]
  grid.points <- cbind(rep(eval.points[,1], ngrid),
  		rep(eval.points[,2], rep(ngrid,ngrid)))
  D <- diff(rbind(hull.points, hull.points[1,]))
  temp  <- D[,1]
  D[,1] <- D[,2]
  D[,2] <- -temp
  C     <- as.vector((hull.points * D) %*% rep(1,2))
  C     <- matrix(rep(C, ngrid^2), nrow = ngrid^2, byrow =T)
  D     <- t(D)
  wy    <- ((grid.points %*% D) >= C)
  wy    <- apply(wy, 1, all)
  wy[wy] <- 1
  wy[!wy] <- NA
  wy    <- matrix(wy, ncol = ngrid)
  }

else {
  w1 <- (w1 > exp(-2))
  w2 <- (w2 > exp(-2))
  wy <- w1 %*% t(w2)
  wy[wy  > 0] <- 1
  wy[wy == 0] <- NA
  }

est <- est * wy

invisible(est)

}



#-----------------Weight matrix for nonparametric regression---------------

sm.weight2 <- function(x, eval.points, h, hmult=1, h.weights=NA, poly.index=1,
				cross = F) {

n  <- nrow(x)
ne <- nrow(eval.points)

if (any(is.na(h.weights))) h.weights <- rep(1,n)

wd1 <-       matrix(rep(eval.points[,1], rep(n, ne)), ncol = n, byrow = T)
wd1 <- wd1 - matrix(rep(x[,1],           ne),         ncol = n, byrow = T)
w   <- matrix(rep(h.weights,             ne),         ncol = n, byrow = T)
w   <- exp(-.5 * (wd1 / (h[1] * hmult * w))^2)
wd2 <-       matrix(rep(eval.points[,2], rep(n, ne)), ncol = n, byrow = T)
wd2 <- wd2 - matrix(rep(x[,2],           ne),         ncol = n, byrow = T)
w   <- w * exp(-.5 * (wd2 / (h[2] * hmult *
	matrix(rep(h.weights, ne), ncol = n, byrow = T)))^2)
if (cross) diag(w) <- 0

if (poly.index==0) {
   den <- w %*% rep(1, n)
   w   <- w / matrix(rep(den,n), ncol=n)
   }

else if (poly.index==1) {
   a11 <-  w              %*% rep(1,n)
   a12 <- (w * wd1      ) %*% rep(1,n)
   a13 <- (w * wd2      ) %*% rep(1,n)
   a22 <- (w * wd1^2    ) %*% rep(1,n)
   a23 <- (w * wd1 * wd2) %*% rep(1,n)
   a33 <- (w * wd2^2    ) %*% rep(1,n)

   d   <- a22 * a33 - a23^2

   b1  <- 1 / (a11 - ((a12*a33 - a13*a23)*a12 + (a13*a22 - a12*a23)*a13)/d)
   b2  <- (a13*a23 - a12*a33) * b1 / d
   b3  <- (a12*a23 - a13*a22) * b1 / d

   wt  <-      matrix(rep(b1, n), ncol = n)
   wt  <- wt + matrix(rep(b2, n), ncol = n) * wd1
   wt  <- wt + matrix(rep(b3, n), ncol = n) * wd2
   w   <- wt * w
   }

w

}

sm.regression.autocor <- function(x=1:n, y, h.first, minh, maxh,
   method="direct", ngrid=15, optimize=F, display="plot", ...)
   {

   GCV<-function(h,x,y,R,sqrt.R,...){
      W<-sm.weight(x,x,h,1,...)
      r <- (y - W %*% as.matrix(y))
      rss<-sum(r^2)
      # ignore autocorrelation
      Trace<-sum(diag(W))
      gcv.0<-rss/(1-Trace/length(x))^2
      # allow for correlation, direct method
      Trace<-sum(diag(W%*%R))
      gcv.r<-rss/(1-Trace/length(x))^2
      # allow for correlation, indirect method
      rw <- backsolve(sqrt.R,r)
      Trace<-sum(diag(W))
      gcv.ri<-sum(rw^2)/(1-Trace/length(x))^2
      c(gcv.0,gcv.r,gcv.ri)
      }
   # end function GCV

   n<-length(y)
   # if(min(diff(x))<=0) stop("x must be in increasing order\n")
   if(length(x)!=n) stop("x and y must have equal length\n")
   if(missing(minh) & missing(x)) minh<-0.5
   if(missing(maxh) & missing(x)) maxh<-10
   w<-sm.weight(x,x,h=h.first,hmult=1)
   ym<-as.vector(w%*%y)
   r<-(y-ym)
   autocov<-rep(0,n)
   for(k in 0:2){
     # estimate first two correlations directly
     u<-r[1:(n-k)]*r[(k+1):n]
     autocov[k+1]<-sum(u)/n
     }
   var<-autocov[1]
   # fit AR(2) model, using first two terms of ACF
   rho1<-autocov[2]/var
   rho2<-autocov[3]/var
   a1 <- rho1*(1-rho2)/(1-rho1^2)
   a2 <- (rho2-rho1^2)/(1-rho1^2)
   type("AR coeff:",c(a1,a2))
   # compute remaining autocorrelations by Yule-Walker equations
   for(k in 3:(n-1)) autocov[k+1]<-a1*autocov[k]+a2*autocov[k-1]
   autocorr<-autocov/var
   R<-diag(n)
   for(k in 1:n) for(j in 1:n) R[k,j]<- autocorr[abs(k-j)+1]
   sqrt.R  <- chol(R)
   hvector <- seq(minh,maxh,length=ngrid)
   min.gcv <- Inf
   h.opt<-0
   result<-matrix(0,ngrid,3,
         dimnames=list(NULL,c("no.cor","direct","indirect")))
   cat("Search for h:")
   for(i in 1:ngrid) {
      h<-hvector[i]
      result[i,]<-GCV(h,x,y,R,sqrt.R)
      cat(" ");cat(i)
           }
   cat("\n")
   # print(cbind(hvector,result))
   if(display=="plot") {
     maxlag<-min(30,n-1)
     acf<-array(autocorr[1:(maxlag+1)],dim=c(maxlag+1,1,1))
     lag<-array(0:maxlag,dim=c(maxlag+1,1,1))
#     acf.plot(list(acf=acf,lag=lag,type="correlation",
#             series="residuals from preliminary smoothing",n.used=n))
        plot(lag, acf, sub="residuals from preliminary smoothing", type="h")
     pause()
     plot(c(hvector[1],hvector[ngrid]),c(min(result),max(result)),
          type="n", xlab="h",
          ylab="Standard and modified GCV criterion",...)
     lines(hvector,result[,method],col=2)
     pause()
     }
   h1<- hvector[order(result[,method])[1]]
   type("Suggested value of h: ",h1)
   sm1<-sm.regression.eval.1d(x,y,h=h1,hmult=1,model="none")
   if(missing(x)) x.name<-"time" else x.name<-deparse(substitute(x))
   if(display=="plot") {
     plot(x,y, xlab=x.name, ylab=deparse(substitute(y)),...)
     lines(sm1$eval.points,sm1$estimate,col=2)
     }
   sm1$aux<-list(h.first=h.first,first.sm=ym,acf=autocorr,raw.residuals=r)
   invisible(sm1)
   }


sm.rm <- function(Time, y, minh=0.1, maxh=2, ngrid=20, optimize=F,
   display="lines", add=F, poly.index=1, display.rice=F, ...){

   rice<-function(h,nSubj,Time,ym,var,r,poly.index=1){
        nTime<-length(Time)
        w<-sm.weight(Time, Time, h, 1, poly.index=poly.index)
        fitted<-w%*%ym
	rss<-sum((ym-fitted)^2)
        Trace<-sum(diag(w%*%r))
        criterion<-rss/nTime-(var/nSubj)*(1-2*Trace/nTime)
	criterion<-sqrt(criterion)
        criterion
       	}

   if(!isMatrix(y)) stop("y must be a matrix")
   nSubj<-dim(y)[1]
   nTime<-dim(y)[2]
   if(missing(Time)) Time <- 1:nTime
   ym<-apply(y,2,mean)
   z <- y-matrix(ym,nrow=nSubj,ncol=nTime,byrow=T)
   autocov<-rep(0,nTime)
   for(k in 0:(nTime-1)){
      u<-z[,1:(nTime-k)]*z[,(k+1):nTime]
      autocov[k+1]<-sum(u)/(nSubj*nTime)
      }
   var<-autocov[1]
   autocorr<-autocov/var
   cat("Autocovariances & autocorrelations:\n")
   print(matrix(cbind(autocov,autocorr),ncol=2,
   dimnames=list(0:(nTime-1),c("auto-cov","auto-corr"))))
   r<-diag(nTime)
   for(k in 1:nTime){
      for(j in 1:nTime) r[k,j]<- autocorr[abs(k-j)+1]
      }
   hvector <-seq(minh,maxh,length=ngrid)
   min.obj<-10^10
   h.opt<-0
   cat("       Rice's criterion:\n")
   cat("     h    indept.   depend.\n")
   result<-matrix(0,ngrid,2,dimnames=list(NULL,c("indept","depend")))
   for(i in 1:ngrid) {
      h<-hvector[i]
      obj.0<-rice(h,nSubj,Time,ym,var,diag(nTime),poly.index)
      obj.r<-rice(h,nSubj,Time,ym,var,r,poly.index)
      result[i,1]<-obj.0
      result[i,2]<-obj.r
      if(obj.r<min.obj){min.obj<-obj.r;  h.opt<-h}
      print(c(h,obj.0,obj.r))
      }
   if(display.rice) {
       plot(c(hvector[1],hvector[ngrid]),c(min(result),max(result)),
         type="n", xlab="h",ylab="sqrt(rice criterion)")
       title(main="Modified Rice criterion for selecting h", sub=
        "\n\n\ndashed line: assume independence, continuous: allow for correlation")
       lines(hvector,result[,1],lty=3)
       lines(hvector,result[,2],lty=1)
       pause()
       }
   if(optimize){
     cat("Search for optimum h using nlminb...\n")
        optimum<-nlminb(start=h.opt,objective=rice,scale=1,lower=0,
         nSubj=nSubj,Time=Time,ym=ym,var=var,r=r,...)
     print(optimum$parameters)
     h.opt<-optimum$parameters
     }
   type("h",h.opt)
   if(display=="se") display1<-"lines" else display1<-display
   sm<-sm.regression(Time,ym,h=h.opt,hmult=1,display=display1, add=add,...)
   if(display=="se"){
      W <- sm.weight(Time,sm$eval.points,h=h.opt,...)
      V <- (var/nSubj) * r
      se <- sqrt(diag(W %*% V %*% t(W)))
      lines(sm$eval.points,sm$estimate+2*se,lty=3)
      lines(sm$eval.points,sm$estimate-2*se,lty=3)
      }
   sm$aux<-list(mean=ym,var=var,autocorr=autocorr,h=h.opt)
   invisible(sm)
   }

sig.trace <- function(expn, hvec, display = "lines")

   {

   expn.char <- paste(deparse(substitute(expn)),collapse="")
   lead.char <- substring(expn.char, 1, nchar(expn.char)-1)
   nvec <- length(hvec)
   pvec <- vector("numeric", length=nvec)

   for (i in 1:nvec) {
      extn.char <- paste(lead.char, ", h = ", as.character(hvec[i]), ")")
      result   <- eval(parse(text = extn.char))
      pvec[i]  <- result$p
      }

   if (!(display == "none")) {
      plot(hvec, pvec, type="l", ylim=c(0,max(pvec)),
      		xlab = "Smoothing parameter, h", ylab = "p-value")
      if (max(pvec) >= 0.05) lines(range(hvec),c(0.05,0.05),lty=2)
      }

   invisible(list(h = hvec, p = pvec))

   }
sm.sphere <- function(lat, long, phi=0, theta=0, kappa=20,
                        panel=F, hidden=F, sphim=F, addpoints=F , ngrid = 32) {

kap <- kappa
invis <- plot2d(lat, long, theta, phi)
sphdraw(theta, phi)

if (!panel) {
        if (hidden) hidplot(invis, theta, phi)
        if (sphim)  sphimage(lat, long, kap, theta, phi, ngrid)
        if (sphim & addpoints)  addplot(lat, long, theta, phi)
        }
else {

items <- c("Set theta and phi",
           "  - increase theta",
           "  - decrease theta",
           "  - increase phi",
           "  - decrease phi",
           "Add hidden points",
           "Add density estimate",
           "  - increase s.p.",
           "  - decrease s.p.",
           "  - add data points",
           "Exit")

ind <- menu(items, graphics=T, title="Sphere")

while (items[ind]!="Exit") {
        if (items[ind]=="Set theta and phi") {
                a <- change(theta, phi)
                theta <- a$theta
                phi <- a$phi
                invis <- plot2d(lat, long, theta, phi)
                sphdraw(theta, phi)
                }
        else if (items[ind]=="  - increase theta") {
                theta <- inctheta(theta, 30)
                invis <- plot2d(lat, long, theta, phi)
                sphdraw(theta, phi)
                }
        else if (items[ind]=="  - decrease theta") {
                theta <- inctheta(theta, -30)
                invis <- plot2d(lat, long, theta, phi)
                sphdraw(theta, phi)
                }
        else if (items[ind]=="  - increase phi") {
                phi <- incphi(phi, 30)
                invis <- plot2d(lat, long, theta, phi)
                sphdraw(theta, phi)
                }
        else if (items[ind]=="  - decrease phi") {
                phi <- incphi(phi, -30)
                invis <- plot2d(lat, long, theta, phi)
                sphdraw(theta, phi)
                }
        else if (items[ind]=="Add hidden points") {
                hidplot(invis, theta, phi)
                }
        else if (items[ind]=="Add density estimate") {
                sphimage(lat, long, kap, theta, phi, ngrid)
                }
        else if (items[ind]=="  - increase s.p.") {
                kap <- kap*2
                sphimage(lat, long, kap, theta, phi, ngrid)
                }
        else if (items[ind]=="  - decrease s.p.") {
                kap <- kap/2
                sphimage(lat, long, kap, theta, phi, ngrid)
                }
        else if (items[ind]=="  - add data points") {
                par(pch="*")
                addplot(lat, long, theta, phi)
                }
        else if (items[ind]=="Add 2nd data set") {
                par(pch="x")
                addplot(lat2, long2, theta, phi)
                }
        ind <- menu(items, graphics=T, title="Sphere")
        }
}

par(pty="m")
invisible(list(theta = theta, phi = phi, kappa = kap))
}

#---------------------------------------------------------------sphdraw

sphdraw<-function(theta, phi)
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

#-------------------------------------------------------------------plot2d

plot2d <- function(d, f, theta, phi)
{
        par(pch="*")
        a <- (f * pi)/180
        b <- (d * pi)/180
        radtheta <- (theta * pi)/180
        radphi <- (phi * pi)/180
        xyzcheck <<- ((cos(a)*cos(b)*sin(radtheta)*cos(radphi)) +
                (sin(b)*sin(radphi)) -
                (sin(a)*cos(b)*cos(radphi)*cos(radtheta)))
        llong <<- a[ xyzcheck >= 0]
        llat <<- b[ xyzcheck >= 0]
        invislong <<- a[ xyzcheck < 0]
        invislat <<- b[ xyzcheck < 0]
        if (length(llat)==0) {
                par(pty="s")
                plot(0,0,type="n", axes=F, xlab="", ylab="", xlim=c(-1.0,1.0),
                        ylim=c(-1.0,1.0))
#               text(-.5,-1.2,labels="Theta=")
#               text(-.2,-1.2,labels=theta)
#               text(0.25,-1.2,labels="Phi=")
#               text(.45,-1.2,labels=phi)
                list(invislong=invislong, invislat=invislat)
                break
                }

        X <<- (cos(llong) * cos(llat) * cos(radtheta)) + (sin(llong) * cos(llat) *
                sin(radtheta))
        Y <<- (sin(llat) * cos(radphi)) + ((cos(llat) * sin(radphi)) * ((sin(llong) *
                cos(radtheta)) - (cos(llong) * sin(radtheta))))
        par(pty = "s")
        plot(X, Y, axes = FALSE, xlab = "", ylab = "", xlim = c(-1.0, 1.0),
                ylim = c(-1.0, 1.0))
#       text(-.5,-1.2,labels="Theta=")
#       text(-.2,-1.2,labels=theta)
#       text(0.25,-1.2,labels="Phi=")
#       text(.45,-1.2,labels=phi)
        list(invislong=invislong, invislat=invislat)
}

#-------------------------------------------------------------------inctheta

inctheta <- function(th, inc)
{

        theta<<-th+inc
        if (theta>=360) theta<<-theta-360
#       cat("Theta =", theta, "\n")
        theta
}

#-------------------------------------------------------------------incphi

incphi <- function(ph, inc)
{
        phi<<- ph+inc
        if (phi>90) phi<<- 90
        if (phi< -90) phi<<- -90
        cat("Phi =", phi,"\n")
        phi
}

#-------------------------------------------------------------------change

change<-function(th, ph)
{
        cat("Theta =",th,"\n")
        cat("Phi =",ph,"\n")
        scan(n=1)
        cat("Change theta to ? \n")
        theta <<- scan(n=1)
        if (theta>=360) theta<<-theta-360
        cat("\n","Change phi to ? \n")
        phi<<- scan(n=1)
        if (phi>90) phi<<- 90
        if (phi< -90) phi<<- -90
        cat("Theta =",theta,"\n")
        cat("Phi =",phi,"\n")
        list(theta=theta, phi=phi)
}

#-------------------------------------------------------------------sphimage

sphimage <- function(latitude, longitude, kap, theta, phi, ngrid = 32)
{
        values <- seq(-1 + 1/ngrid, 1 - 1/ngrid, length = ngrid)
        xgrid  <- rep(values, rep(ngrid, ngrid))
        ygrid  <- rep(values, ngrid)
        dvec   <- rep(0, ngrid**2)
        xlong  <- longitude*pi/180
        xlat   <- latitude*pi/180
        n      <- length(longitude)

        radtheta <- theta*pi/180
        radphi   <- phi*pi/180
        xgrid[xgrid^2 + ygrid^2 >= 1] <- NA
        ygrid[xgrid^2 + ygrid^2 >= 1] <- NA

        za <- -xgrid * sin(radtheta) - ygrid * cos(radtheta) * sin(radphi)
        zb <- cos(radphi) * cos(radtheta) * sqrt(1-xgrid^2-ygrid^2)
        z  <- za + zb
        if ((theta==90) | (theta==270)) x <- -ygrid*sin(radtheta)*sin(radphi) +
                        cos(radphi)*sqrt(1-ygrid^2-z^2)
        else x <- (xgrid + z * sin(radtheta)) / cos(radtheta)
        if (phi==90) y <- sqrt(1-x^2-z^2)
        else if (phi==-90) y <- -sqrt(1-x^2-z^2)
             else y <- (ygrid + (x*sin(radtheta) +
                                z*cos(radtheta))*sin(radphi))/cos(radphi)

        xyzok <- (((x/sqrt(x^2+z^2))*(sqrt(1-y^2))*sin(radtheta)*cos(radphi)) +
                (y*sin(radphi)) -
                   ((-z/sqrt(x^2+z^2))*(sqrt(1-y^2))*cos(radphi)*cos(radtheta)))
        z[xyzok<0] <- (za - zb)[xyzok<0]
        x[xyzok<0] <- ((xgrid + (z*sin(radtheta)))/cos(radtheta))[xyzok<0]
        y[xyzok<0] <- ((ygrid + ((x*sin(radtheta)) +
                        (z*cos(radtheta)))*sin(radphi))/cos(radphi))[xyzok<0]

        xj <- cos(xlong) * cos(xlat)
        yj <- sin(xlat)
        zj <- -sin(xlong) * cos(xlat)
        dvec <- exp(kap*cbind(x,y,z)%*%rbind(xj,yj,zj))%*%rep(1/n,n)
        dvec[is.na(xgrid)] <- 0
        dvec <- dvec / max(dvec)

        fmat <<- matrix(dvec, ngrid, ngrid, byrow = T)

        x <- seq(-1 + 1/ngrid, 1 - 1/ngrid, length = ngrid)
        y <- x
        image(x, y, fmat, add=T)

        angle <- seq(0, pi/2, length = 50)
        xx <- cos(angle)
        yy <- sin(angle)
#	polygon(c(xx,0,1.1,1.1), c(yy,1.1,1.1,0), col=0)
	polygon(c(xx,0,1,1), c(yy,1,1,0), col=0, border=0)
	angle <- seq(pi/2, pi, length = 50)
        xx <- cos(angle)
        yy <- sin(angle)
#	polygon(c(xx,-1.1,-1.1,0),c(yy,0,1.1,1.1), col=0)
	polygon(c(xx,-1,-1,0),c(yy,0,1,1), col=0, border=0)
	angle <- seq(pi, 3*pi/2, length = 50)
        xx <- cos(angle)
        yy <- sin(angle)
#	polygon(c(xx,0,-1.1,-1.1), c(yy,-1.1,-1.1,0), col=0)
	polygon(c(xx,0,-1,-1), c(yy,-1,-1,0), col=0, border=0)
	angle <- seq(3*pi/2, 2*pi, length = 50)
        xx <- cos(angle)
        yy <- sin(angle)
#	polygon(c(xx,1.1,1.1,0), c(yy,0,-1.1,-1.1), col=0)
	polygon(c(xx,1,1,0), c(yy,0,-1,-1), col=0, border=0)
#	polygon(c(-1,-1,1,1), c(-1.4,-1.25,-1.25,-1.4), col=0)
	sphdraw(theta, phi)
}

#-------------------------------------------------------------------addplot

addplot <- function(d, f, theta, phi)
{
        a <- (f * pi)/180
        b <- (d * pi)/180
        radtheta <- (theta * pi)/180
        radphi <- (phi * pi)/180
	xyzcheck <<- ((cos(a)*cos(b)*sin(radtheta)*cos(radphi)) +
(sin(b)*sin(radphi)) - (sin(a)*cos(b)*cos(radphi)*cos(radtheta)))
        llong <<- a[ xyzcheck >= 0]
        llat <<- b[ xyzcheck >= 0]
	if (length(llat)==0) {
			break
		}

        X <<- (cos(llong) * cos(llat) * cos(radtheta)) + (sin(llong) * cos(llat) *
                sin(radtheta))
        Y <<- (sin(llat) * cos(radphi)) + ((cos(llat) * sin(radphi)) * ((sin(llong) *
                cos(radtheta)) - (cos(llong) * sin(radtheta))))
        par(pty = "s")
        points(X, Y)
}

#-------------------------------------------------------------------circle

circle <- function(r)
{
        angle <- seq(0, 7, by = 0.1)
        x <- r * cos(angle)
        y <- r * sin(angle)
	par(lty=1)
        lines(x, y)
}

#-------------------------------------------------------------------hidplot

 hidplot<-function(invis, theta, phi)
{
	invislong <- invis$invislong
	invislat <- invis$invislat
        par(pch = "O")
        a <- (invislong * pi)/180
        b <- (invislat * pi)/180
        radtheta <- (theta * pi)/180
        radphi <- (phi * pi)/180
        if(length(invislat) == 0) {
                points(0, 0, type = "n")
                break
        }
        X <<- (cos(invislong) * cos(invislat) * cos(radtheta)) + (sin(invislong
                ) * cos(invislat) * sin(radtheta))
        Y <<- (sin(invislat) * cos(radphi)) + ((cos(invislat) * sin(radphi)) * (
                (sin(invislong) * cos(radtheta)) - (cos(invislong) * sin(
                radtheta))))
        points(X, Y)
}

#-------------------------------------------------------------------latlines

latlines <- function(beta, theta, phi)
{
	if (beta < (phi-90) | beta >(phi+90)) return()
	par (pch=".")
        radtheta <- (theta * pi)/180
        radbeta <- (beta * pi)/180
        radphi <- (phi * pi)/180
	alpha <- seq(0,(2*pi), by=.05)
	xyzcheck <<- ((cos(alpha)*cos(radbeta)*sin(radtheta)*
cos(radphi)) +(sin(radbeta)*sin(radphi)) - (sin(alpha)*
cos(radbeta)*cos(radphi)*cos(radtheta)))
	alphaplot <- alpha[xyzcheck>=0]
        X <- (cos(alphaplot) * cos(radbeta) * cos(radtheta)) +
(sin(alphaplot) * cos(radbeta) * sin(radtheta))
        Y <- (sin(radbeta) * cos(radphi)) + (((sin(alphaplot) *
cos(radtheta)) - (cos(alphaplot) * sin(radtheta))) * cos(radbeta) *
sin(radphi))
	points(X,Y, pch=".")
}

#-------------------------------------------------------------------latlines.e

latlines.e <- function(beta, theta, phi)
{
	if (beta < (phi-90) | beta >(phi+90)) return()
	par (lty=2)
	par (pch=".")
        radtheta <- (theta * pi)/180
        radbeta <- (beta * pi)/180
        radphi <- (phi * pi)/180
	alpha <- seq(0,(2*pi), by=.005)
	xyzcheck <<- ((cos(alpha)*cos(radbeta)*sin(radtheta)*cos(radphi)) +(sin(radbeta)*sin(radphi)) - (sin(alpha)*cos(radbeta)*cos(radphi)*cos(radtheta)))
	alphaplot <- alpha[xyzcheck>=0]
        X <- (cos(alphaplot) * cos(radbeta) * cos(radtheta)) + (sin(alphaplot) * cos(radbeta) * sin(radtheta))
        Y <- (sin(radbeta) * cos(radphi)) + (((sin(alphaplot) * cos(radtheta)) - (cos(alphaplot) * sin(radtheta))) * cos(radbeta) * sin(radphi))
	points(X,Y, pch=".")
}

#-------------------------------------------------------------------longlines

 longlines <- function(alpha, theta, phi)
{
	par(pch=".")
        radtheta <- (theta * pi)/180
        radalpha <- (alpha * pi)/180
        radphi <- (phi * pi)/180
	beta <- seq(0,(2*pi), by=0.05)
 	xyzcheck <<- ((cos(radalpha) * cos(beta) * sin(radtheta) *
cos(radphi)) + (sin(beta) * sin(radphi)) - (sin(radalpha) * cos(beta) *
cos(radphi) * cos(radtheta)))
	betaplot <- beta[xyzcheck>=0]
     	X <- (cos(radalpha) * cos(betaplot) * cos(radtheta)) +
(sin(radalpha) * cos(betaplot) * sin(radtheta))
      	Y <- (sin(betaplot) * cos(radphi)) + (((sin(radalpha) *
cos(radtheta)) - (cos(radalpha) * sin(radtheta))) * cos(betaplot) *
sin(radphi))
	points(X,Y, pch=".")
}

#-------------------------------------------------------------------longlines.e

 longlines.e <- function(alpha, theta, phi)
{
	par(lty = 2)
	par(pch=".")
        radtheta <- (theta * pi)/180
        radalpha <- (alpha * pi)/180
        radphi <- (phi * pi)/180
	beta <- seq(0,(2*pi), by=0.005)
 	xyzcheck <<- ((cos(radalpha) * cos(beta) * sin(radtheta) * cos(radphi)) + (sin(beta) * sin(radphi)) - (sin(radalpha) * cos(beta) * cos(radphi) * cos(radtheta)))
	betaplot <- beta[xyzcheck>=0]
     	X <- (cos(radalpha) * cos(betaplot) * cos(radtheta)) + (sin(radalpha) * cos(betaplot) * sin(radtheta))
      	Y <- (sin(betaplot) * cos(radphi)) + (((sin(radalpha) * cos(radtheta)) - (cos(radalpha) * sin(radtheta))) * cos(betaplot) * sin(radphi))
	points(X,Y, pch=".")
}

#-------------------------------------------------------------------plot2

plot2 <- function(latitude2, longitude2, theta, phi)

{
        par(pch = "x")
        a <- (longitude2 * pi)/180
        b <- (latitude2 * pi)/180
        radtheta <- (theta * pi)/180
        radphi <- (phi * pi)/180
        xyzcheck <<- ((cos(a) * cos(b) * sin(radtheta) * cos(radphi)) + (sin(b) *
                sin(radphi)) - (sin(a) * cos(b) * cos(radphi) * cos(radtheta)))
        long2 <<- a[xyzcheck >= 0]
        lat2 <<- b[xyzcheck >= 0]
        if(length(lat2) == 0) {
                points(0, 0, type = "n")
                text(0.6, -1.2, labels = "Data set:")
                break
        }
        X <<- (cos(long2) * cos(lat2) * cos(radtheta)) +
(sin(long2) * cos(lat2) * sin(radtheta))
        Y <<- (sin(lat2) * cos(radphi)) + ((cos(lat2)
* sin(radphi)) * ((sin(long2) * cos(radtheta)) - (cos(long2) * sin(radtheta))))
        points(X, Y)
}

#-------------------------------------------------------------------rotate

# rotate <- function (a, b, theta, phi)
# {
# 	theta <<- (theta + a)
#	if (theta>=360) theta <<- (theta-360)
#	if (theta < 0) theta <<- (theta+360)
#	phi <<- (phi + b)
#	if (phi > 90) phi <<- 90
#	if (phi < -90) phi <<- (-90)
#	plot2d(longitude,latitude)
#	sphdraw()
#}

#--------------------------------------------------------------------------
sm.survival <- function(x, y, status, h , hv = 0.05, p = 0.5, status.code = 1,
			eval.points = NA, ngrid = 50, display = "lines",
			xlab = NA, ylab = NA, lty = 1, add = F, ...)

   {

   absent<-function(x) missing(x) | any(is.na(x))
   if (absent(eval.points)) eval.points <- seq(min(x),max(x),length=ngrid)

   if (is.na(xlab))  xlab <- deparse(substitute(x))
   if (is.na(ylab))  ylab <- deparse(substitute(y))

   if (!(display == "none" | add == T)) {
      plot(x, y, type = "n", xlab = xlab, ylab = ylab, ...)
      text(x[status == status.code], y[status == status.code], "x")
      text(x[status != status.code], y[status != status.code], "o")
      }

   n  <- length(x)
   ne <- length(eval.points)

   xr <- x[order(y)]
   statusr <- status[order(y)]
   yr <- sort(y)

   w  <- matrix(rep(eval.points, rep(n, ne)), ncol = n, byrow = T)
   w  <- w - matrix(rep(xr, ne), ncol = n, byrow = T)
   w  <- exp(-.5 * (w/h)^2)

   wf <- t(apply(w,  1, rev))
   wf <- t(apply(wf, 1, cumsum))
   wf <- t(apply(wf, 1, rev))
   w  <- w / wf
   st <- rep(0, n)
   st[statusr == status.code] <- 1
   w  <- 1 - w * matrix(rep(st, ne), ncol = n, byrow= T)

   w  <- w[,st == 1]
   if (ne == 1) w <- matrix(w, ncol = length(w))
   yw <- yr[st == 1]

   w  <- t(apply(w, 1, cumprod))
   w  <- cbind(rep(1, ne), w)
   j  <- -t(apply(w, 1, diff))
   J  <- t(apply(j,  1, cumsum))

   wd <- J - p
   w  <- exp(-.5 * (wd/hv)^2) # * j

   ns <- length(yw)
   s0 <-  w         %*% rep(1,ns)
   s1 <- (w * wd)   %*% rep(1,ns)
   s2 <- (w * wd^2) %*% rep(1,ns)

   w  <-  w * (matrix(rep(s2,ns),ncol=ns) - wd * matrix(rep(s1,ns),ncol=ns))
   w  <-  w / (matrix(rep(s2,ns),ncol=ns) * matrix(rep(s0,ns),ncol=ns)
			- matrix(rep(s1,ns),ncol=ns)^2)
   estimate <- w %*% yw

   if (!(display == "none")) lines(eval.points, estimate, lty = lty)

   invisible(list(estimate = estimate, eval.points = eval.points,
   			h = h, hv = hv))

   }


#  Splus programming tools used by 'sm', but possibly useful elsewhere

provide.data <- function(data, path, describe=TRUE)
{ # load data as data.frame and shows doc file, if it exists
  # assumes that the data are in <data>.dat and documentation in <data>.doc
  # works both on Unix (SunOS) and MS-windows platforms
  name  <- deparse(substitute(data))
  if(missing(path))
     path <- file.path(.sm.home,"data")
  datafile<- file.path(path,paste(name,".dat",sep=""))
  docfile <- file.path(path,paste(name,".doc",sep=""))
  if(!exists(name)){
    if(file.exists(datafile)) {
       cat("Data file being loaded\n")
       assign(name,read.table(datafile,header=T),envir=.GlobalEnv)
       attach(what=data,name=name)  }
    else
       cat("Data file does not exist\n")
      }
  else {
    if(!is.data.frame(data)) cat("object exists, not as a data.frame\n")
    else {
       cat(paste(name,"already loaded\n"))
       attach.frame(data,name=name)}}
  if(describe) file.show(docfile)
  invisible()
} #end function

attach.frame<-function(data, name,...)
{ # attach a data.frame, always in 2nd search position
  if(missing(name)) name<-deparse(substitute(data))
  if(is.data.frame(data)){
    if(!is.na(pos <- match(name, search()))) {
      cat(paste(name,"already attached, re-attached in 2nd position\n"))
      detach(pos=pos)
    }
    cat(paste("attaching",name,"\n",sep=" "))
    attach(what=data, pos=2, name=name,...)}
  else {cat(name);cat(" is not a data.frame\n")}
  invisible()
}


print.graph <- function(file,...){
   dev.print(file=file,...)
   invisible()
}

ask <- function(message="Type in datum")
         eval(parse(prompt=paste(message,": ",sep="")))

type <- function(descr="",x,digits=4)
{  cat(paste(descr," "))
   cat(paste(round(x,digits=digits)))
   cat("\n")
}

pause <- function()
{  cat("Pause. Press <Enter> to continue...")
   readline()
   invisible()
}


sm.ts.pdf <-function(x, h=hnorm(x), lags, maxlag=1,  ask=T, ...)
 {# lag plot of time series, with pdf estimation
   if(missing(lags)) lags<-(1:maxlag) else maxlag<-max(lags)
   if(any(diff(lags))<0) stop("lags must be in increasing order")
   x.name<-deparse(substitute(x))
   x <- as.vector(x)
   n <- length(x)
   marginal<-sm.density(x,ylab="Marginal density", xlab=x.name, ...)
   if(ask) pause()
   for (m in lags){
      x1 <- x[(m+1):n]
      x0 <- x[1:(n-m)]
      biv<-sm.density(cbind(x0,x1), h=rep(h,2),,
             xlab=paste(x.name,"(t-",as.character(m),")",sep=""),
             ylab=paste(x.name,"(t)",sep=""), ...)
      biv$lag <- m
      title(paste("Density of lagged data of ",x.name,
             " (lag=",as.character(m),")",sep=""))
      if(ask & (m<maxlag)) pause()
      }
   invisible(list(marginal=marginal,bivariate=biv))
}


sm.autoregression <-function(x, h=hnorm(x), d=1,maxlag=d, lags,se=F,ask=T, ...)
{# lag plot of time series, with pdf estimation
   x.name<-deparse(substitute(x))
   if(missing(lags)) {
     if (d==1) lags<-(1:maxlag)
     else lags<-cbind(1:(maxlag-1),2:maxlag)
     }
   else
     {if (isMatrix(lags)) d<-2}
   x <- as.vector(x)
   if(d==1) r<-sm.autoregression.1d(x, h, x.name, lags, se=se, ask=ask,...)
   else     r<-sm.autoregression.2d(x, h, x.name, lags, ask=ask)
   invisible(r)
}

sm.autoregression.1d <- function(x,h,x.name,lags,se=F,ask=F,...){
   n <- length(x)
   if(any(diff(lags))<0) stop("lags must be in increasing order")
   x2.name<-paste(x.name,"(t)", sep="")
   xlow<-min(x)-diff(range(x))/20
   xhi<- max(x)+diff(range(x))/20
   lags<-sort(lags)
   for (m in lags){
      x1 <- x[(m+1):n]
      x0 <- x[1:(n-m)]
      r<-sm.regression.eval.1d(x0,x1,h=h,hmult=1,model="none",...)
      x1.name<-paste(x.name,"(t-",as.character(m),")", sep="")
      plot(x0,x1, xlim=c(xlow,xhi), ylim=c(xlow,xhi),
           xlab=x1.name, ylab=x2.name, ...)
      lines(r$eval.points,r$estimate,...)
      if(se) {
        n1 <- length(x0)
        rho1<- cor(x0[-1], x0[-n1])*(1 - 1/n1)
        lines(r$eval.points,r$estimate+2*r$se/sqrt(1-rho1),lty=3)
        lines(r$eval.points,r$estimate-2*r$se/sqrt(1-rho1),lty=3)
        }
      title(paste("Regression of ",x.name, " on past data",sep=""))
      if(ask & (m<lags[length(lags)])) pause()
      }
   invisible(r)
}


sm.autoregression.2d <- function(x, h, x.name, lags, ask=ask,
   ngrid=20, display="none", ...)
  {
   if(dim(lags)[2] != 2) stop("dim(lags)[2] must be 2")
   evpt <- seq(quantile(x,0.1),quantile(x,0.9),length=ngrid)
   n <- length(x)
   nplot<-dim(lags)[1]
   for(ip in 1:nplot){
      m1<-min(lags[ip,])
      m2<-max(lags[ip,])
      x0 <- x[1:(n-m2)]
      x1 <- x[(m2-m1+1):(n-m1)]
      x2 <- x[(m2+1):n]
      r<-sm.regression.eval.2d(cbind(x0,x1),x2, h=c(h,h),hmult=1,
            h.weights=NA,poly.index=1, eval.points=cbind(evpt,evpt), ...)
#      persp(evpt,evpt,r,...)
      persp(evpt,evpt,r,theta = -30, phi = 40, d = 4)
      head<-paste("Regression of ",x.name," on past data (lags: ",
            as.character(m1),", ",as.character(m2),")",sep="")
      title(head)
      if(ask & (ip<nplot)) pause()
      }
   invisible(r)
}
