\name{hcv}
\alias{hcv}
\title{
Cross-validatory choice of smoothing parameter
}
\description{
This function uses the technique of cross-validation to select a smoothing 
parameter suitable for constructing a density estimate or nonparametric regression 
curve in one or two dimensions.
}
\usage{
hcv(x, y=NA, h.weights=NA, ngrid=8, hstart=NA, hend=NA, display="none", add=F, \dots)
}
\arguments{
\item{x}{
a vector, or two-column matrix of data.  If \code{y} is missing these are 
observations to be used in the construction of a density estimate.  If
\code{y} is present, these are the covariate values for a nonparametric regression.
}
\item{y}{
a vector of response values for nonparametric regression.
}
\item{h.weights}{
a vector of weights which multiply the smoothing parameter(s) used in the
kernel function at each observation.
}
\item{ngrid}{
the number of grid points to be used in an initial grid search for the 
value of the smoothing parameter.
}
\item{hstart}{
the smallest value of the grid points to be used in an initial grid search 
for the value of the smoothing parameter.
}
\item{hend}{
the largest value of the grid points to be used in an initial grid search 
for the value of the smoothing parameter.
}
\item{display}{
any character setting other than \code{"none"} will cause the criterion function 
to be plotted over the search grid of smoothing parameters.  The particular
value \code{"log"} will use a log scale for the grid values.
}
\item{add}{
controls whether the plot is added to an existing graph.
}
\item{...}{
additional graphical parameters.
}}
\value{
the value of the smoothing parameter which minimises the cross-validation
criterion over the selected grid.
}
\section{Side Effects}{
If the minimising value is located at the end of the grid of search positions,
or if some values of the cross-validatory criterion cannot be evaluated,
then a warning message is printed.  In these circumstances altering the
values of \code{hstart} and \code{hend} may improve performance.
}
\details{
see Sections 2.4 and 4.5 of the reference below.


The two-dimensional case uses a smoothing parameter derived from a single 
value, scaled by the standard deviation of each component.


This function does not employ a sophisticated algorithm and some adjustment
of the search parameters may be required for different sets of data.
An initial estimate of the value of h which minimises the cross-validatory
criterion is located from a grid search using values which are equally
spaced on a log scale between \code{hstart} and \code{hend}.  A quadratic approximation
is then used to refine this initial estimate.
}
\references{
Bowman, A.W. and Azzalini, A. (1997). Applied Smoothing Techniques for
Data Analysis: the Kernel Approach with S-Plus Illustrations.
Oxford University Press, Oxford.
}
\seealso{
\code{\link{cv}}, \code{\link{hsj}}, \code{\link{hnorm}}
}
\examples{
#  Density estimation


x <- rnorm(50)
par(mfrow=c(1,2))
h.cv <- hcv(x, display="lines", ngrid=32)
sm.density(x, h=hcv(x))
par(mfrow=c(1,1))


#  Nonparametric regression


x <- seq(0, 1, length = 50)
y <- rnorm(50, sin(2 * pi * x), 0.2)
par(mfrow=c(1,2))
h.cv <- hcv(x, y, display="lines", ngrid=32)
sm.regression(x, y, h=hcv(x, y))
par(mfrow=c(1,1))
}
\keyword{nonparametric}
\keyword{smooth}
% Converted by Sd2Rd version 0.3-2.
