\name{nmise}
\alias{nmise}
\title{
mean integrated squared error for density estimation with normal data
}
\description{
This function evaluates the mean integrated squared error of a density
estimate which is constructed from data which follow a normal distribution.
}
\usage{
nmise(sd, n, h)
}
\arguments{
\item{sd}{
the standard deviation of the normal distribution from which the data arise.
}
\item{n}{
the sample size of the data.
}
\item{h}{
the smoothing parameter used to construct the density estimate.
}}
\value{
the mean integrated squared error of the density estimate.
}
\details{
see Section 2.4 of the reference below.
}
\references{
Bowman, A.W. and Azzalini, A. (1997). \emph{Applied Smoothing Techniques for
Data Analysis: the Kernel Approach with S-Plus Illustrations.}
Oxford University Press, Oxford.
}
\seealso{
\code{\link{nise}}
}
\examples{
x  <- rnorm(50)
sd <- sqrt(var(x))
n  <- length(x)
h  <- seq(0.1, 2, length=32)
plot(h, nmise(sd, n, h), type = "l")
}
\keyword{nonparametric}
\keyword{smooth}
% Converted by Sd2Rd version 1.15.
