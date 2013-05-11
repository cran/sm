\name{nise}
\alias{nise}

\title{Integrated squared error between a density estimate and a Normal density}

\description{
This function evaluates the integrated squared error between a density
estimate constructed from a standardised version of the univariate data
\code{y} and a standard normal density function.
}

\usage{nise(y, ...)}

\arguments{
\item{y}{a vector of data.}
\item{...}{further arguments which are to be passed to \code{sm.options}.}
}

\value{
the integrated squared error.
}

\details{
The data \code{y} are first standardised to have sample mean 0 and sample
variance 1.  The integrated squared error between a density estimate
constructed from these standardised data and a standard normal distribution
is then evaluated.

See Section 2.5 of the reference below.
}

\references{
Bowman, A.W. and Azzalini, A. (1997). \emph{Applied Smoothing Techniques for
Data Analysis: the Kernel Approach with S-Plus Illustrations.}
Oxford University Press, Oxford.
}

\seealso{
\code{\link{nmise}}
}

\examples{
x <- rnorm(100)
nise(x)
}

\keyword{nonparametric}
\keyword{smooth}
% Converted by Sd2Rd version 1.15.
