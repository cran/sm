\name{sm.sphere}
\alias{sm.sphere}
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
sm.sphere(lat, long, phi=0, theta=0, kappa=20, panel=F, hidden=F, 
   sphim=F, addpoints=F, ngrid=32)
}
\arguments{
\item{lat}{
a vector giving the latitude component of the data in degrees from the 
equator.
}
\item{long}{
a vector giving the longitude component of the data in degrees east.
}
\item{phi}{
the vertical rotation of the displayed sphere in degrees.  The viewing 
position moves from the equator.
}
\item{theta}{
the horizontal rotation of the displayed sphere.  The viewing position moves
in an easterly direction.  When \code{theta} is set to 0, the origin of the
spherical scales lies in the middle of the right hand side of the plot.
}
\item{kappa}{
the smoothing parameter used to construct the density estimate.  The kernel
function is a Fisher distribution and \code{kappa} is its scale parameter.
Larger values of \code{kappa} will produce smaller amounts of smoothing.
}
\item{panel}{
a logical value which controls the appearance of a panel display.  This
panel may be used to alter a variety of aspects of the estimate and the
plot in an interactive manner.
}
\item{hidden}{
a logical value which, when set to T, will display the points which lie
on the rear side of the displayed sphere.  This argument will be ignored
if \code{sphim} is set to T.
}
\item{sphim}{
a logical value which controls whether a density estimate is constructed
and displayed on the sphere in image form.
}
\item{addpoints}{
a logical value which controls whether the data points are added to the
plot of the density estimate.
}
\item{ngrid}{
the number of points on each side of the regular grid over which the
density estimate is evaluated.
}}
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
Bowman, A.W. and Azzalini, A. (1997). Applied Smoothing Techniques for
Data Analysis: the Kernel Approach with S-Plus Illustrations.
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
sm.sphere(lat, long, sphim=T, kappa=15)
par(mfrow=c(1,1))
}
\keyword{nonparametric}
\keyword{smooth}
% Converted by Sd2Rd version 1.15.
