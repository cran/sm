\name{sm.rm}
\alias{sm.rm}
\title{
Nonparametric analysis of repeated measurements data
}
\description{
This function estimates nonparametrically the mean profile from a matrix
\code{y} which is assumed to contain repeated measurements (i.e. longitudinal
data) from a set of individuals.
}
\usage{
sm.rm(Time, y, minh=0.1, maxh=2, ngrid=20, optimize=F, display="lines",
      add=F, poly.index=1, display.rice=F, \dots)
}
\arguments{
\item{y}{
matrix containing the values of the response variable, with rows associated 
to individuals and columns associated to observation times.
}
\item{Time}{
a vector containing the observation times of the response variable, assumed 
to be the same for all individuals of matrix \code{y}. 
If \code{Time} is not given, this is assumed to be \code{1:ncol(y)}.
}
\item{minh}{
the mimimum value of the interval where the optimal value of the smoothing
parameter is seached according to the modified Rice criterion.
See reference below for details.
}
\item{maxh}{
the maximum value of the above interval.
}
\item{ngrid}{
the number of divisions of the above interval to be considered.
}
\item{optimize}{
Logical value, default is \code{optimize=F}. If \code{optimize=T}, then a full 
optimization is performed after searching the interval \code{(minh,maxh)}
using the optimizer \code{nlminb}.
}
\item{display}{
character value controlling the amount of graphical output of the estimated 
regression curve. It has the same meaning as in \code{sm.regression}. 
Default value is \code{display="lines"}.
}
\item{add}{
logical value, default is \code{add=F}. If \code{add=T} and display is not set
to \code{"none"}, then graphical output added to the existing plot, rather than
starting a new one.
}
\item{poly.index}{
overall degree  of  locally-fitted  polynomial, as used by \code{sm.regression} 
}
\item{display.rice}{
If this set to \code{T} (default is \code{F}), a plot is produced of the curve 
representing the  modified Rice criterion for bandwidth selection. 
See reference below for details.
}
\item{...}{
Optional parameters passed to \code{sm.regression}.
}}
\value{
a list containing the returned value produced by \code{sm.regression} when 
smoothing the mean response value at each given observation time, 
with an extra component \code{$aux} added to the list.
This additional component  is a list itself containing the mean value at each 
observation time, the residual variance of the residuals from the estimated 
regression curve,  the autocorrelation function of the residuals, and the value 
h of the chosen smoothing parameter. 
}
\section{Side Effects}{
if the parameter display is not set to \code{"none"}, a plot of the estimated 
regression curve is produced; 
other aspects are controlled by parameter add and optional parameters (\code{\dots{}}). 
If \code{display.rice=T}, a plot of the modified Rice criterion is shown.
}
\details{
see Section 7.4 of the reference below.
}
\references{
Bowman, A.W. and Azzalini, A. (1997). Applied Smoothing Techniques for
Data Analysis: the Kernel Approach with S-Plus Illustrations.
Oxford University Press, Oxford.
}
\seealso{
\code{\link{sm.regression}}, \code{\link{sm.regression.autocor}}
}
\examples{
provide.data(citrate, describe=FALSE); provide.data(dogs, describe=FALSE)
# assume that citrate and dogs are matrices
a <- sm.rm(y=citrate, display.rice=T)
#
Time <- c(1,3,5,7,9,11,13)
gr1  <- as.matrix(dogs[dogs$Group==1,2:8])
plot(c(1,13), c(3,6),xlab="time", ylab="potassium", type="n") 
sm1  <- sm.rm(Time, gr1, display="se", add=T)
}
\keyword{nonparametric}
\keyword{smooth}
% Converted by Sd2Rd version 0.3-2.
