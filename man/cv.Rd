\name{cv}
\alias{cv}
\title{
Cross-validation criterion for nonparametric density estimation
}
\description{
This function computes a cross-validatory criterion, based on integrated
squared error, for use in selecting a smoothing parameter in nonparametric
density estimation.
}
\usage{
cv(x, h, \dots)
}
\arguments{
\item{x}{
vector, or two-column matrix, of data.
}
\item{h}{
a smoothing parameter.  In the two-dimensional case this is multiplied
by the standard deviation of each component to produce two smoothing
parameters
}
\item{...}{
other optional parameters are passed to the \code{sm.options} function, through
a mechanism which limits their effect only to this call of the function; 
those relevant for this function are the following:
}
\item{h.weights}{
a vector of weights which multiply the smoothing parameter(s) used in the
kernel function at each observation.
}}
\value{
The value of the cross-validatory criterion.
}
\section{Side Effects}{
None
}
\details{
see Section 2.4.3 of the reference below.  The function is called
automatically by \code{hcv} and does not usually need to be called 
independently.
}
\references{
Bowman, A.W. and Azzalini, A. (1997). 
\emph{Applied Smoothing Techniques for Data Analysis: }
\emph{the Kernel Approach with S-Plus Illustrations.}
Oxford University Press, Oxford.
}
\seealso{
\code{\link{hcv}}, \code{\link{hsj}}, \code{\link{hnorm}}, \code{\link{sj}}
}
\examples{
x      <- rnorm(50)
hgrid  <- seq(0.1, 1, length = 10)
cvgrid <- vector("numeric", length = length(hgrid))
for (i in 1:10) cvgrid[i] <- cv(x, hgrid[i])
plot(hgrid, cvgrid, type="l")
}
\keyword{nonparametric}
\keyword{smooth}
% Converted by Sd2Rd version 1.15.
