\name{hnorm}
\alias{hnorm}
\title{
Normal optimal choice of smoothing parameter in density estimation
}
\description{
This functions evaluates the smoothing parameter which is asymptotically
optimal for estimating a density function when the underlying distribution
is Normal.  Data in one, two or three dimensions can be handled.
}
\usage{
hnorm(x)
}
\arguments{
\item{x}{
a vector, or matrix with two or three columns, containing the data.
}}
\value{
the value of the Normal optimal smoothing parameter.
}
\section{Side Effects}{
none
}
\details{
see Section 2.4.2 of the reference below.
}
\references{
Bowman, A.W. and Azzalini, A. (1997). Applied Smoothing Techniques for
Data Analysis: the Kernel Approach with S-Plus Illustrations.
Oxford University Press, Oxford.
}
\seealso{
\code{\link{hcv}}, \code{\link{hsj}}
}
\examples{
x <- rnorm(50)
hnorm(x)
}
\keyword{nonparametric}
\keyword{smooth}
% Converted by Sd2Rd version 0.3-2.
