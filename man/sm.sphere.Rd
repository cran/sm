\name{sm.sphere}
\alias{sm.sphere}
% \alias{latlines}
% \alias{longlines}
% \alias{latlines.e}
% \alias{longlines.e}
% \alias{sphdraw}
% \alias{sphimage}
\title{
Nonparametric density estimation for spherical data.
}
\description{
This function creates a density estimate from data which can be viewed
as lying on the surface of a sphere.  Directional data form a principal
example.  The data are displayed in spherical form and a density estimate
may be superimposed.  The angle of view may be altered.  An interactive 
panel is available to control some features of the estimate and the display.
Only modest amounts of data may be used.  The limit will depend on the
memory available.
}
\usage{
sm.sphere(lat, long, kappa = 20, hidden = FALSE, sphim = FALSE,
          addpoints = FALSE, \dots)
}
\arguments{
  \item{lat}{
    a vector giving the latitude component of the data in degrees from the 
    equator.
  }
  \item{long}{
    a vector giving the longitude component of the data in degrees east.
  }
  \item{kappa}{
    the smoothing parameter used to construct the density estimate.  The kernel
    function is a Fisher distribution and \code{kappa} is its scale parameter.
    Larger values of \code{kappa} will produce smaller amounts of smoothing.
  }
  \item{hidden}{
    a logical value which, when set to \code{TRUE}, will display the points which lie
    on the rear side of the displayed sphere.  This argument will be ignored
    if \code{sphim} is set to \code{TRUE}.
  }
  \item{sphim}{
    a logical value which controls whether a density estimate is constructed
    and displayed on the sphere in image form.
  }
  \item{addpoints}{
    a logical value which controls whether the data points are added to the
    plot of the density estimate.
  }
  \item{\dots}{
    arguments for \code{\link{sm.options}}.
  }
}
\value{
a list containing the value of the smoothing parameter and the rotation 
angles of the displayed plot.
}
\section{Side Effects}{
none.
}
\details{
see Section 1.5 of the reference below.
}
\references{
Bowman, A.W. and Azzalini, A. (1997). \emph{Applied Smoothing Techniques for
Data Analysis: the Kernel Approach with S-Plus Illustrations.}
Oxford University Press, Oxford.
}
\seealso{
\code{\link{sm.density}}
}
\examples{
lat  <- rnorm(50, 10, 15)
long <- c(rnorm(25, 300, 15), rnorm(25, 240, 15))
par(mfrow=c(1,2))
sm.sphere(lat, long)
sm.sphere(lat, long, sphim=TRUE, kappa=15)
par(mfrow=c(1,1))
}
\keyword{nonparametric}
\keyword{smooth}
% Converted by Sd2Rd version 1.15.
