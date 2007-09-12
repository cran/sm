\name{sm.variogram}
\alias{sm.variogram}
\title{
A test of spatial independence based on smoothing an empirical variogram.
}
\description{
  This function constructs an empirical variogram, using the robust form
  of construction based on square-root absolute value differences of the
  data.  Nonparametric regression is used to assess whether the mean
  value of the empirical variogram changes with distance.  A hypothesis
  test and a graphical reference band are both available, under the
  assumption that the errors in the data are approximately normally
  distributed.
}
\usage{
sm.variogram(x, y, h, \dots)
}
\arguments{
\item{x}{
  a vector or two-column matrix of spatial location values.
}
\item{y}{
  a vector of responses observed at the spatial locations.
}
\item{h}{
a smoothing parameter to be used on the distance scale.  A normal kernel
function is used and \code{h} is its standard deviation.  However, if
this argument is omitted \code{h} will be selected by an approximate
degrees of freedom criterion, controlled by the \code{df} parameter.
See \code{sm.options} for details.
}
\item{\dots}{
  other optional parameters are passed to the \code{sm.options}
  function, through a mechanism which limits their effect only to this
  call of the function; those relevant for this function are 
  \code{add},
  \code{eval.points}, 
  \code{ngrid}, 
  \code{se}, 
  \code{xlab}, 
  \code{ylab}, 
  \code{xlim}, 
  \code{ylim}, 
  \code{lty};
  see the documentation of  \code{\link{sm.options}} for their
  description.  See the details section below for a discussion of
  the \code{display} and \code{se} parameters in this setting.
}}  
\value{
  a list containing the raw (\code{sqrtdiff}, \code{distance}) and binned 
  (\code{sqrtdiff.mean}, \code{distance.mean}) differences and distances,
  the frequencies of the bins (\code{weights}), 
  the values of the estimate (\code{estimate}) at the evaluation points 
  (\code{eval.points}) the value of the smoothing parameter (\code{h}) used,
  an indicator (\code{ibin}) of the bin in which the distance between
  each pair of observations was placed and the indices (\code{ipair}) of
  the original observations used to construct each pair.  
  If a test of spatial independence has been performed, the p-value (\code{p}) 
  is returned.  If \code{se} was set to \code{TRUE} then the standard errors
  (\code{se}) of the binned values, under the assumption of spatial independence, 
  are returned.  If a reference band has been requested, the standard error 
  (\code{se.band}) of the difference between the smooth estimate and the mean
  of all the data points is also returned.
}
\section{Side Effects}{
  a plot on the current graphical device is produced, unless the option 
  \code{display="none"} is set.  
}
\details{
  The reference below describes the statistical methods used in the function.

  The \code{display} argument has a special meaning for
  this function.  Its default value is \code{"binned"}, which plots the
  binned version of the empirical variogram.  As usual, the value \code{"none"}
  will suppress the graphical display.  Any other value will lead to a plot 
  of the individual differences between all observations.  This will lead
  to a very large number of plotted points, unless the dataset is small.
  
  Setting \code{se = TRUE} produces standard errors for the binned values of the 
  empirical variogram.  They are constructed to account for the correlation which
  arises from the differencing involved in the variogram.  However, the calculations
  assume spatial independence and they therefore simply give an indication of the 
  variability under this assumption.  
}
\references{
  Diblasi, A. and Bowman, A.W. (2001).
  On the use of the variogram for checking independence in a 
  Gaussian spatial process.
  \emph{Biometrics}, 57, 211-218.
}
\seealso{
\code{\link{sm.regression}}, \code{\link{sm.options}}
}
\examples{
provide.data(coalash)
Position <- cbind(East, North)
sm.options(df = 6)

par(mfrow=c(2,2))
sm.variogram(Position, Percent, band = FALSE)
sm.variogram(Position, Percent)
sm.variogram(East,     Percent)
sm.variogram(North,    Percent)
par(mfrow=c(1,1))
}
\keyword{smooth}
\keyword{regression}
\keyword{spatial}
