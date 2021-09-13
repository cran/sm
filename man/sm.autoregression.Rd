\name{sm.autoregression}
\alias{sm.autoregression}
\title{
Nonparametric estimation of the autoregression function
}
\description{
This function estimates nonparametrically the autoregression function
(conditional mean given the past values) of a time series \code{x},
assumed to be stationary.
}
\usage{
sm.autoregression(x, h = hnorm(x), d = 1, maxlag = d, lags,
                  se = FALSE, ask = TRUE)
}
\arguments{
\item{x}{
vector containing the time series values.
}
\item{h}{
the bandwidth used for kernel smoothing.
}
\item{d}{
number of past observations used for conditioning; it must be 1
(default value) or 2.
}
\item{maxlag}{
maximum of the lagged values to be considered (default value is \code{d}).
}
\item{lags}{
if \code{d==1}, this is a vector containing the lags considered for conditioning;
if \code{d==2}, this is a matrix with two columns, whose rows contains pair of
values considered for conditioning.
}
\item{se}{
if \code{se==T}, pointwise confidence bands are computed of approximate level 95\%.
}
\item{ask}{
if \code{ask==TRUE}, the program pauses after each plot until <Enter> is pressed.
}
}
\value{
a list with the outcome of the final estimation (corresponding to
the last value or pairs of values of lags), as returned by \code{sm.regression}.
}
\section{Side Effects}{
graphical output is produced on the current device.
}
\details{
see Section 7.3 of the reference below.
}
\references{
Bowman, A.W. and Azzalini, A. (1997). 
\emph{Applied Smoothing Techniques for Data Analysis: }
\emph{the Kernel Approach with S-Plus Illustrations.}
Oxford University Press, Oxford.
}
\seealso{
\code{\link{sm.regression}}, \code{\link{sm.ts.pdf}}
}
\examples{
sm.autoregression(log(lynx), maxlag=3, se=TRUE)
sm.autoregression(log(lynx), lags=cbind(2:3,4:5))
}
\keyword{nonparametric}
\keyword{smooth}
\keyword{ts}
% Converted by Sd2Rd version 1.15.
