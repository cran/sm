\name{sm.density.compare}
\alias{sm.density.compare}
\title{
Comparison of univariate density estimates
}
\description{
This function allows a set of univariate density estimates to be
compared, both graphically and formally in a bootstrap hypothesis test of 
equality.
}
\usage{
sm.density.compare(x, group, h=NA, model="none", test=TRUE,
                   nboot=100, monitor=TRUE, \dots)
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
the default value is \code{"none"} which restricts comparison to plotting only.
The alternative value \code{"equal"} can produce a bootstrap hypothesis test
of equality and the display of an appropriate reference band.
}
\item{test}{
a logical flag controlling the production of a bootstrap test of equality.
}
\item{band}{
a logical flag controlling the production of a reference band for equality.  
A band will be produced only in the case of two groups.
}
\item{nboot}{
the number of bootstrap simulations.
}
\item{monitor}{
a logical flag controlling the printing of the iteration numbers during the
bootstrap test.
}
\item{\dots}{
additional \code{sm.options} or graphical parameters.
}}
\value{
When \code{model} is set to \code{"none"}, nothing is returned.  When
\code{"model"} is set to \code{"equal"}, a list containing the smoothing
parameter and the p-value of the test is returned.  When \code{band}
takes the value \code{TRUE}, and there are only two groups to compare, the list
contains in addition the upper and lower end-points of the reference
band for equality.
}
\section{Side Effects}{
none.
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
\code{\link{sm.density}}, \code{\link{sm.ancova}}
}
\examples{
y <- rnorm(100)
g <- rep(1:2, rep(50,2))
sm.density.compare(y, g, model="equal")
}
\keyword{nonparametric}
\keyword{smooth}
% Converted by Sd2Rd version 1.15.
