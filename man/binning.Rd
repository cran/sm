\name{binning}
\alias{binning}
\title{
Construct frequency table from raw data
}
\description{
Given a vector or a matrix \code{x}, this function constructs a frequency table
associated to appropriate intervals covering the range of x.
}
\usage{
binning(x, breaks, nbins)
}
\arguments{
\item{x}{
a vector or a matrix with either one or two columns. 
If \code{x} is a one-dimentional matrix, this is equivalent to a vector.
}
\item{breaks}{
either a vector or a matrix with two columns, assigning the division
points of the axis, or the axes in the matrix case.
If \code{breaks} is not given, it is computed by dividing the range of \code{x}
into \code{nbins} intervals for each of the axes.
}
\item{nbins}{
the number of intervals on the \code{x} axis (in the vector case), 
or a vector of two elements with the number of intervals on each
axes of \code{x} (in the marix case).
If \code{nbins} is not given, a value is computed as \code{round(log(length(x),2)+1)}
or using a similar expression in the matrix case.
}}
\value{
in the vector case, this is a list containing the vector \code{midpoints} 
of the interval midpoints and the frequecies \code{freq} associated to them;
in the matrix case, the returned value is a list with the following 
elements: a two-dimensional matrix \code{x} with the coordinates of the
midpoints of the two-dimensional bins excluding those with 0 frequecies, 
its associated matrix \code{x.freq} of frequencies, the coodinateds of the 
\code{midpoints}, the division points, and the observed frequencies \code{freq.table}
in full tabular form.
}
\details{
This function is primarity intended for use in connection with \code{sm.density},
to estimate noparametrically a density function, when the number of data 
points is high.  To avoid lengthy computations and use of very large matrices, 
the data are tabulated with the use of \code{binning}, and the outcome is
passed to \code{sm.density  } which computes the estimated density curve,
using methods described in Chapter 1 of the reference below.
}
\references{
Bowman, A.W. and Azzalini, A. (1997). Applied Smoothing Techniques for
Data Analysis: the Kernel Approach with S-Plus Illustrations.
Oxford University Press, Oxford.
}
\seealso{
\code{\link{sm.density}}, \code{\link{cut}},\code{\link{table}}
}
\examples{
# example of 1-d use
x<-rnorm(1000)
xb<-binning(x)
sm.density(xb$x,h=hnorm(x),weights=xb$freq)
# example of 2-d use
x<-rnorm(1000)
x<-cbind(x,x+rnorm(1000))
xb<-binning(x)
plot(x)
sm.density(xb$x, h=hnorm(x), weights=xb$x.freq, display="slice", add=T)
}
% Converted by Sd2Rd version 0.3-2.
