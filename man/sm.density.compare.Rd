\name{sm.density.compare}
\alias{sm.density.compare}
\title{
Comparison of univariate density estimates
}
\description{
This function allows a set of univariate density estimates to be
compared, both graphically and formally in a permutation test of 
equality.
}
\usage{
sm.density.compare(x, group, h, model = "none",  \dots)
}
\arguments{
\item{x}{
a vector of data.
}
\item{group}{
a vector of group labels.  If this is not already a factor it will be converted to a factor.
}
\item{h}{
the smoothing parameter to be used in the construction of each density
estimate.  Notice that the same smoothing parameter is used for each group.
If this value is omitted, the mean of the normal optimal values for the
different groups is used.
}
\item{model}{
the default value is \code{"none"} which restricts comparison to
plotting only. The alternative value \code{"equal"} can produce a
bootstrap hypothesis test of equality and the display of an appropriate
reference band.
}
\item{...}{
other optional parameters are passed to the \code{sm.options} function,
through a mechanism which limits their effect only to this call of the
function. Relevant parameters for this function are:
\code{method}, \code{df}, \code{band}, \code{test}, \code{nboot},
plus those controlling graphical display (unless \code{display = "none"}
is set) such as \code{col}, \code{col.band}, \code{lty} and \code{lwd};
see the documentation of  \code{\link{sm.options}} for their description.
The parameter \code{nboot} controls the number of permutations used in the
permutation test.
}}
\value{
A list is returned containing:
\item{estimate}{a matrix whose rows contain the estimates for each group.}
\item{eval.points}{the grid of common evaluation points for the estimates.}
\item{h}{the common smoothing parameter used in the construction of the estimates.}
\item{levels}{the levels of the group factor.}
\item{col, lty, lwd}{plotting details which can be useful in constructing a legend for the plot; see the examples below.}

When \code{"model"} is set to \code{"equal"}, the output list also contains the p-value (\code{p}) of the test.  

When \code{band = TRUE}, and there are only two groups to compare, the output list also contains the upper (\code{upper}) and lower (\code{lower}) edges of the reference band for equality.
}

\section{Side Effects}{
 a plot on the current graphical device is produced, unless
  \code{display="none"}.
}
\details{
For a general description of  the methods involved, see Section 6.2 of the reference below.

The colours and linetypes of the density estimates are set by \code{col} and \code{lty} which can be passed as additional arguments.  By default these are set to \code{1 + 1:ngroup}, where \code{ngroup} is the number of groups represented in the \code{group} variable.
}
\references{
Bowman, A.W. and Azzalini, A. (1997). \emph{Applied Smoothing Techniques for
Data Analysis: the Kernel Approach with S-Plus Illustrations.}
Oxford University Press, Oxford.
}
\seealso{
\code{\link{sm.density}}, \code{\link{sm.ancova}}, \code{\link{sm.options}}
}
\examples{
y <- rnorm(100)
g <- rep(1:2, rep(50,2))
sm.density.compare(y, g)

comp <- sm.density.compare(y, g, model = "equal")
legend("topleft", comp$levels, col = comp$col, lty = comp$lty, lwd = comp$lwd)
}
\keyword{nonparametric}
\keyword{smooth}
% Converted by Sd2Rd version 1.15.
