\name{binning}
\alias{binning}
\title{
Construct frequency table from raw data
}
\description{
Given a vector or a matrix \code{x}, this function constructs a frequency table
associated to appropriate intervals covering the range of \code{x}.
}
\usage{
binning(x, y, breaks, nbins)
}
\arguments{
\item{x, y}{
a vector or a matrix with either one or two columns. 
If \code{x} is a one-dimentional matrix, this is equivalent to a vector.
}
\item{breaks}{
either a vector or a matrix with two columns (depending on the dimension of \code{x}),
assigning the division points of the axis, or the axes in the matrix case.
It must not include \code{Inf},\code{-Inf} or \code{NA}s, and it must span the whole range of 
the \code{x} points.
If \code{breaks} is not given, it is computed by dividing the range of \code{x}
into \code{nbins} intervals for each of the axes.
}
\item{nbins}{
the number of intervals on the \code{x} axis (in the vector case), 
or a vector of two elements with the number of intervals on each
axes of \code{x} (in the matrix case).
If \code{nbins} is not given, a value is computed as \code{round(log(length(x))/log(2)+1)}
or using a similar expression in the matrix case.
}}
\value{
In the vector case, a list is returned containing the following elements:
a vector \code{x} of the midpoints of the bins excluding those with 0 frequecies, 
its associated matrix \code{x.freq} of frequencies, the coodinateds of the 
\code{midpoints}, the division points, and the complete vector of observed 
frequencies \code{freq.table} (including the 0 frequencies), and the vector
\code{breaks} of division points.
In the matrix case, the returned value is a list with the following 
elements: a two-dimensional matrix \code{x} with the coordinates of the
midpoints of the two-dimensional bins excluding those with 0 frequecies, 
its associated matrix \code{x.freq} of frequencies, the coordinates of the 
\code{midpoints}, the matrix \code{breaks} of division points, and the observed 
frequencies \code{freq.table} in full tabular form.
}
\details{
This function is called automatically (under the default settings)
by some of the functions of the \code{sm} library when the sample size is
large, to allow handling of datasets of essentially unlimited size.
Specifically, it is used by \code{sm.density}, \code{sm.regression}, \code{sm.ancova},
\code{sm.binomial} and \code{sm.poisson}.
}
\references{
Bowman, A.W. and Azzalini, A. (1997). 
\emph{Applied Smoothing Techniques for Data Analysis:
the Kernel Approach with S-Plus Illustrations.}
Oxford University Press, Oxford.
}
\seealso{
\code{\link{sm}}, \code{\link{sm.density}}, \code{\link{sm.regression}}, \code{\link{sm.binomial}}, \code{\link{sm.poisson}}, \code{\link{cut}}, \code{\link{table}}
}
\examples{
# example of 1-d use
x  <- rnorm(1000)
xb <- binning(x)
xb <- binning(x, breaks=seq(-4,4,by=0.5))
# example of 2-d use
x <- rnorm(1000)
y <- 2*x + 0.5*rnorm(1000)
x <- cbind(x, y)
xb<- binning(x, nbins=12)
}
\keyword{nonparametric}
% Converted by Sd2Rd version 1.15.
