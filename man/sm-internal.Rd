\name{sm-internal}
\alias{addplot}
\alias{attach.frame}
\alias{box1.}
\alias{box2.}
\alias{britmap}
\alias{change}
\alias{circle}
\alias{coord2.}
\alias{hidplot}
\alias{incphi}
\alias{inctheta}
\alias{normdens.band}
\alias{np.contour.plot.3d.}
\alias{p.quad.moment}
\alias{plot.regression}
\alias{plot2}
\alias{plot2d}
\alias{print.graph}
\alias{replace.na}
\alias{sm.glm}
\alias{sm.imageplot}
\alias{sm.persplot}
\alias{sm.sigma}
\alias{sm.sigweight}
\alias{sm.sliceplot}
\alias{sm.weight}
\alias{sm.weight2}
\alias{smplot.density}
\alias{type}
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
attach.frame(data, name, ...)
box1.(theta = pi/6, phi = pi/6, sc, col = par("col"), axes.lim = sc)
box2.(theta = pi/6, phi = pi/6, sc, labels = c("", "", ""),
    col = par("col"), cex = 9/10, axes.lim = sc)
britmap()
change(th, ph)
circle(r)
coord2.(x, y, z, theta, phi, sc)
hidplot(invis, theta, phi)
incphi(ph, inc)
inctheta(th, inc)
isInteger(x)
isMatrix(x)
normdens.band(x, h, weights = rep(1, length(x)), options = list())
np.contour.plot.3d.(coord, data = matrix(), shadow = TRUE, gridsize = 20,
    numpts = 3, xmin = NA, xmax = NA, ymin = NA, ymax = NA, zmin = NA,
    zmax = NA, xlab = NA, ylab = NA, zlab = NA, theta = pi/4,
    phi = pi/4, colour = FALSE, title.colour = 3, label.colour = 3,
    axes.colour = 6, plot.colour = "blue", shadow.colour = "cyan",
    cex = NA)
p.quad.moment(A, Sigma, tobs, ndevs)
plot.regression(x, y, design.mat, h, r, model, weights, rawdata = list(),
    options = list(), ...)
plot2(latitude2, longitude2, theta, phi)
plot2d(d, f, theta, phi)
print.graph(file, ...)
replace.na(List, comp, value)
sm.glm(x, y, family, h, eval.points, start, offset)
sm.imageplot(x, y, h, weights, rawdata, options = list())
sm.persplot(x, y, h = hnorm(cbind(x, y), weights), weights, rawdata = list(),
    options = opt)
sm.sigma(x, y, rawdata, weights = rep(1, length(y)), diff.ord = 2)
sm.sigweight(x, weights = rep(1, length(x)), ...)
sm.sliceplot(x, y, h, weights, rawdata = list(), options = list())
sm.weight(x, eval.points, h, cross = FALSE, weights = rep(1, n), options)
sm.weight2(x, eval.points, h, cross = FALSE, weights = rep(1, nrow(x)),
    options = list())
smplot.density(x, h, weights = rep(1, length(x)), rawdata = list(x = x),
    options = list())
type(descr = "", x, digits = 4)
wmean(x, w)
wvar(x, w)
}
\details{
  These are not to be called by the user.
}
\keyword{internal}
