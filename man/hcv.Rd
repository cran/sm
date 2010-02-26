\name{hcv}
\alias{hcv}
\title{
Cross-validatory choice of smoothing parameter
}
\description{
This function uses the technique of cross-validation to select a smoothing 
parameter suitable for constructing a density estimate or nonparametric 
regression curve in one or two dimensions.
}
\usage{
hcv(x, y = NA, hstart = NA, hend = NA, \dots) 
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
\item{hstart}{
the smallest value of the grid points to be used in an initial grid search 
for the value of the smoothing parameter.
}
\item{hend}{
the largest value of the grid points to be used in an initial grid search 
for the value of the smoothing parameter.
}
\item{\dots}{
   other optional parameters are passed to the \code{sm.options} function,
through a mechanism which limits their effect only to this call of the
function. Those specifically relevant for this function are the following:
\code{h.weights}, \code{ngrid}, \code{display}, \code{add};
see the documentation of  \code{\link{sm.options}} for their description.
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
See Sections 2.4 and 4.5 of the reference below.

The two-dimensional case uses a smoothing parameter derived from a single 
value, scaled by the standard deviation of each component.

This function does not employ a sophisticated algorithm and some
adjustment of the search parameters may be required for different sets
of data.  An initial estimate of the value of h which minimises the
cross-validatory criterion is located from a grid search using values
which are equally spaced on a log scale between \code{hstart} and
\code{hend}.  A quadratic approximation is then used to refine this
initial estimate.
}
\note{As from version 2.1 of the package, a similar effect can be
  obtained with the new function \code{h.select}, via
  \code{h.select(x, method="cv")}. Users are encouraged to adopt
  this route, since \code{hcv} might be not accessible directly
  in future releases of the package.  When the
  sample size is large \code{hcv} uses the raw data while
  \code{h.select(x, method="cv")} uses binning.  The latter is 
  likely to produce a more stable choice for \code{h}.
}
\references{
Bowman, A.W. and Azzalini, A. (1997).
\emph{Applied Smoothing Techniques for Data Analysis:} 
\emph{the Kernel Approach with S-Plus Illustrations.}
Oxford University Press, Oxford.
}
\seealso{
   \code{\link{h.select}}, \code{\link{hsj}}, \code{\link{hnorm}}
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
% Converted by Sd2Rd version 1.15.
