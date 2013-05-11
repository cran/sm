\name{sm-internal}
\alias{addplot}
\alias{aicc.crit.reg}
% \alias{ask}
% \alias{attach.frame}
\alias{britmap}
\alias{change}
\alias{circle}
\alias{cv}
\alias{cv.crit.dens}
\alias{cv.crit.reg}
\alias{df.crit.reg}
\alias{hidplot}
\alias{incphi}
\alias{inctheta}
\alias{normdens.band}
\alias{np.contour.plot.3d.}
\alias{p.quad.moment}
\alias{plot2}
\alias{plot2d}
\alias{replace.na}
\alias{sj}
\alias{sm.density.1d}
\alias{sm.density.2d}
\alias{sm.density.3d}
\alias{sm.density.eval.1d}
\alias{sm.density.eval.2d}
\alias{sm.density.positive.1d}
\alias{sm.density.positive.2d}
\alias{sm.density.positive.grid}
\alias{sm.glm}
\alias{sm.imageplot}
\alias{sm.persplot}
\alias{sm.rglplot}
\alias{sm.regression.1d}
\alias{sm.regression.2d}
\alias{sm.regression.eval.1d}
\alias{sm.regression.eval.2d}
\alias{sm.regression.test}
\alias{sm.sigweight}
\alias{sm.sliceplot}
\alias{sm.weight}
\alias{sm.weight2}
\alias{smplot.density}
\alias{smplot.regression}
\alias{wmean}
\alias{wvar}
\alias{isInteger}
\alias{isMatrix}
\title{Internal sm functions}
\description{
 Internal \code{sm} functions
}
\usage{
addplot(d, f, theta, phi)
britmap()
change(th, ph)
circle(r)
cv(x, h, ...)
hidplot(invis, theta, phi)
incphi(ph, inc)
inctheta(th, inc)
isInteger(x)
isMatrix(x)
normdens.band(x, h, weights = rep(1, length(x)), options = list())
p.quad.moment(A, Sigma, tobs, ndevs)
smplot.regression(x, y, design.mat, h, r, model, weights, rawdata = list(),
    options = list(), ...)
plot2(latitude2, longitude2, theta, phi)
plot2d(d, f, theta, phi)
replace.na(List, comp, value)
sj(x, h)
sm.check.data(x, y = NA, weights = NA, group = NA, ...)
sm.density.1d(x, h = hnorm(x, weights), model = "none", weights,
   rawdata = list(x = x), options = list())
sm.density.2d(X, h = hnorm(X, weights), weights = rep(1, length(x)),
   rawdata = list(), options = list())
sm.density.3d(x, h = hnorm(x, weights),  weights = rep(1, length(x)), 
   rawdata = list(), options = list())
sm.density.eval.1d(x, h, weights = rep(1, n), options = list())
sm.density.eval.2d(x, y, h, xnew, ynew, eval.type = "points",
   weights = rep(1, n), options = list())
sm.density.positive.1d(x, h, weights, options = list())
sm.density.positive.2d(X, h = c(hnorm(log(X[, 1] + delta[1]), weights),
   hnorm(log(X[,2] + delta[2]), weights)), eval.type = "points",
   weights = rep(1, nrow(X)), options = list())
sm.density.positive.grid(X, h = c(hnorm(log(X[, 1] + delta[1])),
   hnorm(log(X[, 2] + delta[2]))), weights=NA, options=list())
sm.glm(x, y, family, h, eval.points, start, offset, options=list())
sm.imageplot(x, y, h, weights, rawdata, options = list())
sm.persplot(x, y, h = hnorm(cbind(x, y), weights), weights, rawdata = list(),
    options = list())
sm.regression.1d(x, y, h, design.mat = NA, model = "none",
    weights = rep(1, length(x)), rawdata, options = list())
sm.regression.2d(x, y, h, model = "none", weights = rep(1, length(y)), rawdata,
    options = list())
sm.regression.eval.1d(x, y, design.mat, h, model = "none",
    weights = rep(1, length(x)), rawdata, options = list())
sm.regression.eval.2d (x, y, h, model, eval.points, hull = TRUE,
    weights, options = list())
sm.regression.test(x, y, design.mat = NA, h, model = "no.effect",
    weights = rep(1,length(y)), rawdata, options = list())
sm.sigweight(x, weights)
sm.sliceplot(x, y, h, weights, rawdata = list(), options = list())
sm.weight(x, eval.points, h, cross = FALSE, weights = rep(1, length(x)), options)
sm.weight2(x, eval.points, h, cross = FALSE, weights = rep(1, nrow(x)),
    options = list())
smplot.density(x, h, weights = rep(1, length(x)), rawdata = list(x = x),
    options = list())
wmean(x, w)
wvar(x, w)
}
\details{
  These are not to be called by the user.
}
\keyword{internal}
