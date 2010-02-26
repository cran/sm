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
a vector of group labels.
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
function. Those specifically relevant for this function are the following:
\code{method}, \code{df}, \code{band}, \code{test}, \code{nboot},
plus those controlling graphical display (unless \code{display="none"}
is set);
see the documentation of  \code{\link{sm.options}} for their description.
The parameter \code{nboot} controls teh number of permutations used in the
permutation test.
}}
\value{
When \code{model} is set to \code{"none"}, nothing is returned.  When
\code{"model"} is set to \code{"equal"}, a list containing the smoothing
parameter and the p-value of the test is returned.  When \code{band}
takes the value \code{TRUE}, and there are only two groups to compare,
the list contains in addition the upper and lower end-points of the reference
band for equality.
}
\section{Side Effects}{
 a plot on the current graphical device is produced, unless
  \code{display="none"}.
}
\details{
see Section 6.2 of the reference below.
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
sm.density.compare(y, g, model="equal")
}
\keyword{nonparametric}
\keyword{smooth}
% Converted by Sd2Rd version 1.15.
