\name{sm.discontinuity}

\alias{sm.discontinuity}

\title{The detection of discontinuities in a regression curve or surface.
}

\description{
  This function uses a comparison of left and right handed nonparametric 
  regression curves to assess the evidence for the presence of one or
  more discontinuities in a regression curve or surface.  A hypothesis
  test is carried out, under the assumption that the errors in the data 
  are approximately normally distributed.  A graphical indication of the 
  locations where the evidence for a discontinuity is strongest is also
  available.
}

\usage{
sm.discontinuity(x, y, h, hd, \dots)
}
\arguments{
\item{x}{
  a vector or two-column matrix of covariate values.
}
\item{y}{
  a vector of responses observed at the covariate locations.
}
\item{h}{
a smoothing parameter to be used in the construction of the nonparametric
regression estimates.  A normal kernel
function is used and \code{h} is its standard deviation(s).  However, if
this argument is omitted \code{h} will be selected by an approximate
degrees of freedom criterion, controlled by the \code{df} parameter.
See \code{sm.options} for details.
}
\item{hd}{
a smoothing parameter to be used in smoothing the differences of the
left and right sided nonparametric regression estimates.  A normal kernel
function is used and \code{hd} is its standard deviation(s).  However, if
this argument is omitted \code{hd} will be set to \code{h * sqrt(0.25)}, and
\code{h} reset to \code{h * sqrt(0.75)}, when \code{x} is a vector
When \code{x} is a matrix, \code{hd} will be set to \code{h * sqrt(0.5)}
and \code{h} will be reset to the same value.
}
\item{\dots}{
  other optional parameters are passed to the \code{sm.options}
  function, through a mechanism which limits their effect only to this
  call of the function; those relevant for this function are 
  \code{add},
  \code{eval.points}, 
  \code{ngrid}, 
  \code{se}, 
  \code{band},
  \code{xlab}, 
  \code{ylab}, 
  \code{xlim}, 
  \code{ylim}, 
  \code{lty},
  \code{col};
  see the documentation of  \code{\link{sm.options}} for their
  description.
}}  
\value{
  a list containing the following items
  \item{p}{the p-value for the test of the null hypothesis that no
           discontinuities are present.}
  \item{sigma}{the estimated standard deviation of the errors.}
  \item{eval.points}{the evaluation points of the nonparametric
           regression estimates.  When \code{x} is a matrix, 
           \code{eval.points} is also a matrix whose columns
           define the evaluation grid of each margin of the 
           evaluation rectangle.}
  \item{st.diff}{a vector or matrix of standardised differences
           between the left and right sided estimators at the 
           evaluation points.}
  \item{diffmat}{when \code{x} is a vector, this contains the
           locations and standardised differences where the latter
           are greater than 2.5.}
  \item{angle}{when \code{x} is a matrix, this contains the estimated 
           angles at which the standardised differences were constructed.}
  \item{h}{the principal smoothing parameter.}
  \item{hd}{the smoothing parameter used for double-smoothing (see the
           reference below).}
}
\section{Side Effects}{
  a plot on the current graphical device is produced, unless the option 
  \code{display="none"} is set.  
}
\details{
  The reference below describes the statistical methods used in the function.
  There are minor differences in some computational details of the implementation.
  
  Currently duplicated rows of \code{x} cause a difficulty in the two covariate case.  Duplicated rows should be removed.
}
\references{
   Bowman, A.W., Pope, A. and Ismail, B. (2006).
   Detecting discontinuities in nonparametric regression curves and
   surfaces.
   \emph{Statistics & Computing}, 16, 377--390.
}
\seealso{
\code{\link{sm.regression}}, \code{\link{sm.options}}
}
\examples{
par(mfrow = c(3, 2))

with(nile, {
   sm.discontinuity(Year, Volume, hd = 0)
   sm.discontinuity(Year, Volume)

   ind <- (Year > 1898)
   plot(Year, Volume)
   h <- h.select(Year, Volume)
   sm.regression(Year[!ind], Volume[!ind], h, add = TRUE)
   sm.regression(Year[ ind], Volume[ ind], h, add = TRUE)

   hvec <- 1:15
   p <- numeric(0)
   for (h in hvec) {
      result <- sm.discontinuity(Year, Volume, h,
                          display = "none", verbose = 0)
      p <- c(p, result$p)
   }
   plot(hvec, p, type = "l", ylim = c(0, max(p)), xlab = "h")
   lines(range(hvec), c(0.05, 0.05), lty = 2)
})

with(trawl, {
   Position  <- cbind(Longitude, Latitude)
   ind <- (Longitude < 143.8)
   # Remove a repeated point which causes difficulty with sm.discontinuity
   ind[54] <- FALSE
   sm.regression(Position[ind,], Score1[ind], theta = 35, phi = 30)
   sm.discontinuity(Position[ind,], Score1[ind], col = "blue")
})
par(mfrow = c(1, 1))
	
#  The following example takes longer to run.
#  Alternative values for nside are 32 and 64.
#  Alternative values of yjump are 1 and 0.5.
# nside  <- 16
# yjump  <- 2
# x1     <- seq(0, 1, length = nside)
# x2     <- seq(0, 1, length = nside)
# x      <- expand.grid(x1, x2)
# x      <- cbind(x1 = x[, 1], x2 = x[, 2])
# y      <- rnorm(nside * nside)
# ind    <- (sqrt((x[, 1] - 0.5)^2 + (x[, 2] - 0.5)^2) <= 0.25)
# y[ind] <- y[ind] + yjump
# image(x1, x2, matrix(y, ncol = nside))
# sm.discontinuity(x, y, df = 20, add = TRUE)
}
\keyword{smooth}
\keyword{regression}
