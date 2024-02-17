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
sm.rm(Time, y, minh = 0.1, maxh = 2, optimize = FALSE,
      rice.display = FALSE, \dots)
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
the minimum value of the interval where the optimal value of the smoothing
parameter is searched according to the modified Rice criterion.
See reference below for details.
}
\item{maxh}{
the maximum value of the above interval.
}
\item{optimize}{
  Logical value, default is \code{optimize=FALSE}. If
  \code{optimize=TRUE}, then a full 
  optimization is performed after searching the interval \code{(minh,maxh)}
  using the optimizer \code{optim}.
}
\item{rice.display}{
  If this set to \code{TRUE} (default is \code{FALSE}), a plot is
  produced of the curve 
  representing the modified Rice criterion for bandwidth selection. 
  See reference below for details.
}
\item{\dots}{
  other optional parameters are passed to the \code{sm.options}
  function, through a mechanism which limits their effect only to this
  call of the function; those relevant for this function are the following:
  \describe{
    \item{add}{
      logical value, default is \code{add=FALSE}. If \code{add=TRUE} and
      display is not set to \code{"none"}, then graphical output added
      to the existing plot, rather than starting a new one.
    }
    \item{display}{
      character value controlling the amount of graphical output of the estimated 
      regression curve. It has the same meaning as in \code{sm.regression}. 
      Default value is \code{display="lines"}.
    }
    \item{ngrid}{
      the number of divisions of the above interval to be considered. 
      Default: \code{ngrid=20}.
    }
    \item{poly.index}{
      overall degree  of  locally-fitted  polynomial, as used by
      \code{sm.regression}. Default: \code{ngrid=1}.
    }
  }
}
}
\value{
a list containing the returned value produced by \code{sm.regression} when 
smoothing the mean response value at each given observation time, 
with an extra component \code{$aux} added to the list.
This additional component  is a list itself containing the mean value at each 
observation time, the residual variance of the residuals from the estimated 
regression curve,  the autocorrelation function of the residuals, and 
the value \code{h} of the chosen smoothing parameter. 
}
\section{Side Effects}{
if the parameter display is not set to \code{"none"}, a plot of the estimated 
regression curve is produced; 
other aspects are controlled by the optional parameters (\code{\dots}). 
If \code{rice.display=TRUE}, a plot of the modified Rice criterion is shown.
}
\details{
see Section 7.4 of the reference below.
}
\references{
Bowman, A.W. and Azzalini, A. (1997). 
\emph{Applied Smoothing Techniques for Data Analysis: }
\emph{the Kernel Approach with S-Plus Illustrations.}
Oxford University Press, Oxford.
}
\seealso{
\code{\link{sm.regression}}, \code{\link{sm.regression.autocor}}, \code{\link{optim}}
}
\examples{
sm.rm(y=as.matrix(citrate), display.rice=TRUE)
#
with(dogs, {
   Time <- seq(1,13,by=2)
   gr1  <- as.matrix(dogs[dogs$Group==1,2:8])
   plot(c(1,13), c(3,6),xlab="time", ylab="potassium", type="n") 
   sm1  <- sm.rm(Time, gr1, display="se", add=TRUE)
})
}
\keyword{nonparametric}
\keyword{smooth}
% Converted by Sd2Rd version 1.15.
