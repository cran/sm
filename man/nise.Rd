\name{nise}
\alias{nise}
\title{
integrated squared error between a density estimate and a Normal density
}
\description{
This function evaluates the integrated squared error between a density
estimate constructed from a standardised version of the univariate data
\code{y} and a standard normal density function.
}
\usage{
nise(y, hmult=1)
}
\arguments{
\item{y}{
a vector of data.
}
\item{hmult}{
a factor which can be used to multiply the normal optimal smoothing parameter
before construction of the density estimate.
}}
\value{
the integrated squared error.
}
\section{Side Effects}{
none.
}
\details{
the data \code{y} are first standardised to have sample mean 0 and sample
variance 1.  The integrated squared error between a density estimate
constructed from these standardised data and a standard normal distribution
is then evaluated.


see Section 2.5 of the reference below.
}
\references{
Bowman, A.W. and Azzalini, A. (1997). Applied Smoothing Techniques for
Data Analysis: the Kernel Approach with S-Plus Illustrations.
Oxford University Press, Oxford.
}
\seealso{
\code{\link{nmise}}
}
\examples{
x <- rnorm(x)
nise(x)
}
\keyword{nonparametric}
\keyword{smooth}
% Converted by Sd2Rd version 0.3-2.
