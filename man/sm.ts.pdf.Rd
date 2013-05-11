\name{sm.ts.pdf}
\alias{sm.ts.pdf}
\title{
Nonparametric density estimation of stationary time series data
}
\description{
This function estimates the density function of a time series \code{x},
assumed to be stationary. The univariate marginal density is estimated
in all cases; bivariate densities of pairs of lagged values are estimated
depending on the parameter \code{lags}.
}
\usage{
sm.ts.pdf(x, h = hnorm(x), lags, maxlag = 1, ask = TRUE)
}
\arguments{
\item{x}{
a vector containing a time series
}
\item{h}{
bandwidth
}
\item{lags}{
  for each value, \code{k} say, in the vector \code{lags} a density
  estimate is produced
  of the joint distribution of the pair \code{(x(t-k),x(t))}.
}
\item{maxlag}{
  if \code{lags} is not given, it is assigned the value \code{1:maxlag}
  (default=1).
}
\item{ask}{
  if \code{ask=TRUE}, the program pauses after each plot, until <Enter>
  is pressed.
}
}
\value{
a list of two elements, containing the outcome of the estimation of 
the marginal density and the last bivariate density, as produced by 
\code{\link{sm.density}}.
}
\section{Side Effects}{
plots are produced on the current graphical device.
}
\details{
see Section 7.2 of the reference below.
}
\references{
Bowman, A.W. and Azzalini, A. (1997). \emph{Applied Smoothing Techniques for
Data Analysis: the Kernel Approach with S-Plus Illustrations.}
Oxford University Press, Oxford.
}
\seealso{
\code{\link{sm.density}}, \code{\link{sm.autoregression}}
}
\examples{
with(geyser, {
   sm.ts.pdf(geyser$duration, lags=1:2)
})
}
\keyword{nonparametric}
\keyword{smooth}
\keyword{ts}
% Converted by Sd2Rd version 1.15.
