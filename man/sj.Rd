\name{sj}
\alias{sj}
\title{
Sheather-Jones criterion for nonparametric density estimation
}
\description{
This function computes a criterion associated with the Sheather-Jones
plug-in method of selecting a smoothing parameter in nonparametric
density estimation.  The selected smoothing parameter is identified
by the point at which this criterion takes the value 0.
}
\usage{
sj(x, h)
}
\arguments{
\item{x}{
a vector of data.
}
\item{h}{
a value of smoothing parameter.
}}
\value{
The value of the Sheather-Jones criterion.
}
\section{Side Effects}{
none.
}
\details{
see Section 2.4.4 of the reference below.  The function is called
automatically by \code{hsj} and does not usually need to be called 
independently.
}
\references{
Bowman, A.W. and Azzalini, A. (1997). Applied Smoothing Techniques for
Data Analysis: the Kernel Approach with S-Plus Illustrations.
Oxford University Press, Oxford.
}
\seealso{
\code{\link{cv}}, \code{\link{hcv}}, \code{\link{hsj}}, \code{\link{hnorm}}
}
\examples{
x      <- rnorm(50)
hgrid  <- seq(0.1, 1, length = 10)
sjgrid <- vector("numeric", length = length(hgrid))
for (i in 1:10) sjgrid[i] <- sj(x, hgrid[i])
plot(hgrid, sjgrid, type="l")
abline(0, 0, lty=2)
}
\keyword{nonparametric}
\keyword{smooth}
% Converted by Sd2Rd version 0.3-2.
